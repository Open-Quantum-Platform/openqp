/* Consumed by tests/test_mp2_dispatch.py; fail one real Fortran allocation.
 * dyld rewrites malloc calls from loaded images to this function. The shim's
 * own call remains bound to libSystem malloc, per dyld interpose semantics. */
#include <stddef.h>
#include <stdlib.h>
#include <stdatomic.h>
#include <errno.h>
static _Atomic size_t target_bytes;
static _Atomic int countdown;
static _Atomic int failures;
void mp2_test_arm(size_t bytes, int occurrence) {
    atomic_store(&failures, 0);
    atomic_store(&countdown, occurrence);
    atomic_store(&target_bytes, bytes);
}
int mp2_test_failures(void) { return atomic_load(&failures); }
static void *mp2_test_malloc(size_t bytes) {
    if (bytes == atomic_load(&target_bytes) && atomic_load(&countdown) > 0) {
        if (atomic_fetch_sub(&countdown, 1) == 1) {
            atomic_fetch_add(&failures, 1);
            errno = ENOMEM;
            return NULL;
        }
    }
    return malloc(bytes);
}
__attribute__((used)) static struct { const void *replacement; const void *original; }
mp2_malloc_interpose __attribute__((section("__DATA,__interpose"))) = {
    (const void *)mp2_test_malloc, (const void *)malloc
};
