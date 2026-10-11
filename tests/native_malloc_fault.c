/* Consumed by test_nevpt2_subspaces.py: glibc-only malloc fault injection. */
#define _GNU_SOURCE
#include <stddef.h>
#include <link.h>
#include <errno.h>
#include <stdint.h>
extern void *__libc_malloc(size_t);
extern void *__libc_calloc(size_t, size_t);
/* This shim is loaded at process startup through LD_PRELOAD. Allocator
 * callbacks must use static TLS: general-dynamic TLS calls __tls_get_addr,
 * which can allocate while glibc grows the DTV for a late-loaded module and
 * recursively re-enter malloc. Initial-exec accesses the startup TLS directly.
 */
#define STARTUP_TLS _Thread_local __attribute__((tls_model("initial-exec")))
static STARTUP_TLS size_t target_size;
static STARTUP_TLS int remaining = -1;
static STARTUP_TLS int failures;
static STARTUP_TLS uintptr_t owner_start, owner_end;

/* Resolve the executable mapping while unarmed, never in an allocator callback.
 * BLAS may allocate the same number of bytes as a Fortran scratch array. Only
 * allocations directly requested by the native owner count toward the ordinal.
 */
static int find_owner(struct dl_phdr_info *obj, size_t size, void *symbol) {
    uintptr_t address = (uintptr_t)symbol;
    (void)size;
    for (size_t i = 0; i < obj->dlpi_phnum; ++i) {
        const ElfW(Phdr) *ph = &obj->dlpi_phdr[i];
        uintptr_t start = obj->dlpi_addr + ph->p_vaddr;
        uintptr_t end = start + ph->p_memsz;
        if (ph->p_type == PT_LOAD && (ph->p_flags & PF_X) &&
                address >= start && address < end) {
            owner_start = start; owner_end = end;
            return 1;
        }
    }
    return 0;
}
int oqp_test_set_malloc_owner(void *symbol) {
    remaining = -1;
    owner_start = owner_end = 0;
    return dl_iterate_phdr(find_owner, symbol);
}
void oqp_test_arm_malloc(size_t size, int ordinal) {
    target_size = size; remaining = ordinal; failures = 0;
}
int oqp_test_malloc_failures(void) { return failures; }
static int fail_allocation(size_t size, void *caller) {
    uintptr_t address = (uintptr_t)caller;
    if (remaining >= 0 && address >= owner_start && address < owner_end &&
            size == target_size && remaining-- == 0) {
        ++failures; errno = ENOMEM; return 1;
    }
    return 0;
}
void *malloc(size_t size) {
    if (fail_allocation(size, __builtin_return_address(0))) return NULL;
    return __libc_malloc(size);
}
void *calloc(size_t count, size_t size) {
    if (size && count <= SIZE_MAX / size && fail_allocation(count * size, __builtin_return_address(0)))
        return NULL;
    return __libc_calloc(count, size);
}
