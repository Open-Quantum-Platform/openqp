#define _GNU_SOURCE
/* Consumed by tests/test_mp2_dispatch.py; fail one real Fortran allocation.
 * dyld rewrites malloc calls from loaded images to this function. The shim's
 * own call remains bound to libSystem malloc, per dyld interpose semantics. */
#include <stddef.h>
#include <stdlib.h>
#include <stdatomic.h>
#include <errno.h>
#include <string.h>
static _Atomic size_t target_bytes;
static _Atomic int countdown;
static _Atomic size_t start_bytes;
static _Atomic int failures;
static _Atomic int select_n5;
static _Atomic int selector_calls;
/* The production selector reads this existing cap once to choose a path,
 * and again only after failure to decide whether direct-J is affordable.
 * Control those two queries, on every rank, without changing nbf > 120 or
 * any numerical operation. No allocator is called from this getenv hook. */
void mp2_test_select_n5(void) {
    atomic_store(&selector_calls, 0);
    atomic_store(&select_n5, 1);
}
int mp2_test_selector_calls(void) { return atomic_load(&selector_calls); }
#ifdef __APPLE__
static char *mp2_test_getenv(const char *name)
#else
char *getenv(const char *name)
#endif
{
    if (atomic_load(&select_n5) && strcmp(name, "OQP_MP2_MAX_JBUILDS") == 0) {
        return atomic_fetch_add(&selector_calls, 1) == 0 ? "0" : "4000";
    }
#ifdef __APPLE__
    return getenv(name);
#else
    /* Use the original environment directly; dlsym would itself allocate. */
    extern char **environ;
    size_t length = strlen(name);
    for (char **entry = environ; *entry; ++entry) {
        if (strncmp(*entry, name, length) == 0 && (*entry)[length] == '=')
            return *entry + length + 1;
    }
    return NULL;
#endif
}
void mp2_test_arm(size_t bytes, int occurrence, size_t start) {
    atomic_store(&failures, 0);
    atomic_store(&start_bytes, start);
    atomic_store(&countdown, occurrence);
    atomic_store(&target_bytes, bytes);
}
int mp2_test_failures(void) { return atomic_load(&failures); }
#ifdef __APPLE__
static void *mp2_test_malloc(size_t bytes)
#else
void *malloc(size_t bytes)
#endif
{
    if (atomic_load(&start_bytes) == bytes) atomic_store(&start_bytes, 0);
    if (atomic_load(&start_bytes) == 0 && bytes == atomic_load(&target_bytes) && atomic_load(&countdown) > 0) {
        if (atomic_fetch_sub(&countdown, 1) == 1) {
            atomic_fetch_add(&failures, 1);
            errno = ENOMEM;
            return NULL;
        }
    }
#ifdef __APPLE__
    return malloc(bytes);
#else
    /* Linux/glibc: the un-interposed libc allocator. */
    extern void *__libc_malloc(size_t);
    return __libc_malloc(bytes);
#endif
}
#ifdef __APPLE__
__attribute__((used)) static struct { const void *replacement; const void *original; }
mp2_interpose[] __attribute__((section("__DATA,__interpose"))) = {
    { (const void *)mp2_test_malloc, (const void *)malloc },
    { (const void *)mp2_test_getenv, (const void *)getenv }
};

#endif
