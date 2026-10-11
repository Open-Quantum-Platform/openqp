/* Exercise the TagArray C ABI in the installed liboqp, not inline Record code. */
#define _POSIX_C_SOURCE 200809L
#include <dlfcn.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

/* Public TagArray 1.0 C ABI (tagarray.h and tagarray/RecordInfo.h). */
typedef struct {
    int32_t type_id, itemsize;
    int64_t count, ndims;
    int64_t *dims;
    void *data;
    char *desc;
} RecordInfo;

static void *symbol(void *library, const char *name) {
    void *address = dlsym(library, name);
    if (!address) {
        fprintf(stderr, "%s: %s\n", name, dlerror());
        exit(2);
    }
    return address;
}

#define CHECK(condition) do { \
    if (!(condition)) { \
        fprintf(stderr, "TagArray check failed at line %d: %s\n", __LINE__, #condition); \
        exit(3); \
    } \
} while (0)

/* Also exercise the production cofactor caller and its debug records. */
static void cofactor_lifetime(void *library) {
    void (*persistent)(int64_t, int64_t, const int64_t *, const double *,
                        double *, double *, int *, int *) =
        symbol(library, "occupied_cofactors_persistent_selftest");
    const int64_t sizes[4] = {6, 7, 3, 6};
    double inputs[7 * 7 * 4] = {0.0};
    double actual[7 * 7 * 4], previous[7 * 7 * 4];
    for (int step = 0; step < 4; ++step) {
        for (int i = 0; i < sizes[step]; ++i) {
            inputs[step * 49 + i * 7 + i] = 1.0;
            if (i + 1 < sizes[step])
                inputs[step * 49 + (i + 1) * 7 + i] = 0.125;
        }
    }
    CHECK(setenv("NAC_DUMP_MINORS", "1", 1) == 0);
    for (int repeat = 0; repeat < 3; ++repeat) {
        int checks[6] = {-1, -1, -1, -1, -1, -1}, status = -1;
        persistent(7, 4, sizes, inputs, actual, previous, checks, &status);
        printf("Cofactor status=%d checks=%d,%d,%d,%d,%d,%d\n", status,
               checks[0], checks[1], checks[2], checks[3], checks[4], checks[5]);
        CHECK(fflush(stdout) == 0);
        CHECK(status == 0);
        CHECK(checks[0] == 1 && checks[1] == 4 && checks[2] == 16);
        CHECK(checks[3] == 4 && checks[4] == 4 && checks[5] == 2);
        for (int step = 0; step < 4; ++step) {
            for (int j = 0; j < 7; ++j) {
                for (int i = 0; i < 7; ++i) {
                    double expected = 0.0;
                    if (i < sizes[step] && j < sizes[step] && i >= j) {
                        expected = 1.0;
                        for (int k = j; k < i; ++k) expected *= -0.125;
                    }
                    int offset = step * 49 + j * 7 + i;
                    CHECK(actual[offset] == expected);
                    CHECK(previous[offset] == expected);
                }
            }
        }
        const int64_t invalid[4] = {6, 1, 7, 3};
        persistent(7, 4, invalid, inputs, actual, previous, checks, &status);
        CHECK(status == 1 && checks[0] == 0 && checks[5] == 2);
        CHECK(unsetenv("NAC_DUMP_MINORS") == 0);
        persistent(7, 4, sizes, inputs, actual, previous, checks, &status);
        CHECK(status == 4 && checks[0] == 1 && checks[1] == 1);
        CHECK(checks[2] == 0 && checks[5] == 2);
        CHECK(setenv("NAC_DUMP_MINORS", "1", 1) == 0);
    }
    puts("Cofactor same-owner lifetime regression passed (3 x 4 sizes)");
}

int main(int argc, char **argv) {
    CHECK(argc == 2);
    void *library = dlopen(argv[1], RTLD_NOW | RTLD_LOCAL);
    if (!library) {
        fprintf(stderr, "dlopen: %s\n", dlerror());
        return 2;
    }
    void *(*new_container)(const char *) = symbol(library, "TA_Container_new");
    void (*delete_container)(void *) = symbol(library, "TA_Container_delete");
    int32_t (*create)(void *, const char *, int32_t, int32_t, const int64_t *,
                      const uint8_t *, const char *) = symbol(library, "TA_Container_create");
    int32_t (*append)(void *, const char *, int32_t, const int64_t *,
                      const uint8_t *) = symbol(library, "TA_Container_append");
    RecordInfo (*get)(void *, const char *) = symbol(library, "TA_Container_get");
    void (*erase)(void *, const char *) = symbol(library, "TA_Container_erase");
    void (*clear)(void *) = symbol(library, "TA_Container_clear");

    const int64_t sizes[] = {1, 7, 64, 257, 4096, 19, 2};
    double input[4096];
    for (int i = 0; i < 4096; ++i) input[i] = i + 0.25;
    for (int repeat = 0; repeat < 3; ++repeat) {
        void *container = new_container("linked lifetime regression");
        CHECK(container);
        clear(container);  /* empty cleanup is valid */
        for (unsigned s = 0; s < sizeof(sizes) / sizeof(sizes[0]); ++s) {
            int64_t n = sizes[s];
            /* TA_TYPE_REAL64 = 10; one owned record and one live sibling. */
            CHECK(create(container, "record", 10, 1, &n,
                         (const uint8_t *)input, "") == 0);
            CHECK(create(container, "sibling", 10, 1, &n, NULL, "") == 0);
            RecordInfo record = get(container, "record");
            CHECK(record.count == n && record.itemsize == sizeof(double));
            CHECK((uintptr_t)record.data % 64 == 0);
            CHECK(((double *)record.data)[n - 1] == input[n - 1]);
            CHECK(append(container, "record", 1, &n, (const uint8_t *)input) == 0);
            record = get(container, "record");  /* previous view was invalidated */
            CHECK(record.count == 2 * n);
            CHECK(((double *)record.data)[2 * n - 1] == input[n - 1]);
            erase(container, "record");
            RecordInfo sibling = get(container, "sibling");
            CHECK(sibling.count == n && ((double *)sibling.data)[n - 1] == 0.0);
            clear(container);
            clear(container);  /* repeated cleanup must not double-free */
        }
        /* Leave a record for the container destructor to release. */
        int64_t n = 13;
        CHECK(create(container, "owned_at_exit", 10, 1, &n, NULL, "") == 0);
        delete_container(container);
    }
    cofactor_lifetime(library);
    CHECK(dlclose(library) == 0);
    puts("TagArray linked lifetime regression passed (3 x 7 sizes)");
    return 0;
}
