/* Exercise the TagArray C ABI in the installed liboqp, not inline Record code. */
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
    CHECK(dlclose(library) == 0);
    puts("TagArray linked lifetime regression passed (3 x 7 sizes)");
    return 0;
}
