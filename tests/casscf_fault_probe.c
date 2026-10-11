/* Linked by source/CMakeLists.txt; used by casscf_driver/casscf_anhess and
 * tests/test_casscf_status.py to enter production Hessian/final-CI failures. */
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include "../include/oqp.h"

#ifdef _MSC_VER
#define CAS_TLS __declspec(thread)
#else
#define CAS_TLS _Thread_local
#endif
static CAS_TLS int32_t enabled, selected, hits;

int32_t oqp_casscf_test_fault_begin(int32_t site)
{
    const char *gate = getenv("OQP_CASSCF_TEST_FAULTS");
    const char *test = getenv("PYTEST_CURRENT_TEST");
    if (!gate || strcmp(gate, "1") || !test ||
        !strstr(test, "test_casscf_status.py::"))
        return 1;
    if (enabled)
        return 2;
    if (site != 1 && site != 2 && site != 3)
        return 3;
    enabled = 1;
    selected = site;
    hits = 0;
    return 0;
}

int32_t oqp_casscf_test_fault_end(int32_t *count)
{
    if (!count)
        return 1;
    *count = hits;
    enabled = selected = hits = 0;
    return 0;
}

int32_t oqp_casscf_test_fault_request(int32_t site)
{
    if (!enabled || !selected || site != selected)
        return 0;
    selected = 0; /* Consume once, before Python retry or native recovery. */
    hits++;
    return 1;
}
