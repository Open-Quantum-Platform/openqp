/* Consumed by test_nevpt2_subspaces.py: glibc-only malloc fault injection. */
#include <stddef.h>
#include <errno.h>
#include <stdint.h>
extern void *__libc_malloc(size_t);
extern void *__libc_calloc(size_t, size_t);
static _Thread_local size_t target_size;
static _Thread_local int remaining = -1;
static _Thread_local int failures;
void oqp_test_arm_malloc(size_t size, int ordinal) {
    target_size = size; remaining = ordinal; failures = 0;
}
int oqp_test_malloc_failures(void) { return failures; }
static int fail_allocation(size_t size) {
    if (remaining >= 0 && size == target_size && remaining-- == 0) {
        ++failures; errno = ENOMEM; return 1;
    }
    return 0;
}
void *malloc(size_t size) {
    if (fail_allocation(size)) return NULL;
    return __libc_malloc(size);
}
void *calloc(size_t count, size_t size) {
    if (size && count <= SIZE_MAX / size && fail_allocation(count * size))
        return NULL;
    return __libc_calloc(count, size);
}
