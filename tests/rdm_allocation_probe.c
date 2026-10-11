/* Linked by source/CMakeLists.txt; called by fci_sigma_strings.F90 and
 * tests/test_rdm_fallback.py. No compiler allocator or private symbol ABI. */
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include "../include/oqp.h"

struct probe_state {
    int32_t enabled;
    int32_t fail_ordinal;
    struct oqp_rdm_test_alloc_stats stats;
};

/* Windows ifx/icx builds use the Microsoft DLL TLS spelling; other supported
 * C compilers provide C11 TLS. RDM allocations occur on the calling thread
 * before its OpenMP workers start, including in OpenMP-OFF builds. */
#ifdef _MSC_VER
static __declspec(thread) struct probe_state probe;
#else
static _Thread_local struct probe_state probe;
#endif

int32_t oqp_rdm_test_alloc_begin(int32_t ordinal)
{
    const char *enable = getenv("OQP_RDM_TEST_ALLOCATIONS");
    const char *test = getenv("PYTEST_CURRENT_TEST");
    if (!enable || strcmp(enable, "1") || !test ||
        (!strstr(test, "test_rdm_fallback.py::") &&
         !strstr(test, "rdm_posix_allocator.py::test_rdm_never_changes_foreign_blas_policy")))
        return 1;
    if (probe.enabled)
        return 2;
    if (ordinal < 0 || ordinal > 18)
        return 3;
    memset(&probe, 0, sizeof(probe));
    probe.enabled = 1;
    probe.fail_ordinal = ordinal;
    return 0;
}

int32_t oqp_rdm_test_alloc_end(struct oqp_rdm_test_alloc_stats *stats)
{
    if (!stats)
        return 1;
    *stats = probe.stats;
    memset(&probe, 0, sizeof(probe));
    return 0;
}

int32_t oqp_rdm_test_alloc_active(void)
{
    return probe.enabled;
}

/* Called on the controller's owning OS thread after the OpenMP barrier.
 * Workers never access the controller's TLS. */
void oqp_rdm_test_contraction_record(int32_t team)
{
    if (probe.enabled && team > probe.stats.max_contraction_team)
        probe.stats.max_contraction_team = team;
}

int32_t oqp_rdm_test_alloc_request(int64_t n1, int64_t n2, int64_t n3,
                                 int32_t element_bytes)
{
    uint64_t bytes;
    if (!probe.enabled)
        return 0;
    if (n1 < 0 || n2 < 0 || n3 < 0 || element_bytes <= 0)
        return 1;
    bytes = (uint64_t)n1;
    if ((n2 && bytes > (uint64_t)INT64_MAX / (uint64_t)n2))
        return 1;
    bytes *= (uint64_t)n2;
    if ((n3 && bytes > (uint64_t)INT64_MAX / (uint64_t)n3))
        return 1;
    bytes *= (uint64_t)n3;
    if (bytes > (uint64_t)INT64_MAX / (uint64_t)element_bytes)
        return 1;
    bytes *= (uint64_t)element_bytes;
    if (bytes > SIZE_MAX || probe.stats.count == INT32_MAX ||
        bytes > (uint64_t)(INT64_MAX - probe.stats.total_bytes))
        return 1;
    probe.stats.count++;
    probe.stats.total_bytes += (int64_t)bytes;
    if ((int64_t)bytes > probe.stats.largest_bytes)
        probe.stats.largest_bytes = (int64_t)bytes;
    if (probe.fail_ordinal && probe.stats.count == probe.fail_ordinal) {
        probe.fail_ordinal = 0; /* Consume once; recovery needs no rearming. */
        return 1;
    }
    return 0;
}
