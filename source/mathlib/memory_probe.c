/*
 * How much memory this process may actually use.
 *
 * Modules that allocate a large array up front (conventional coupled
 * cluster, for example) want to refuse before the allocator does, with a
 * message that says what was needed and what was there.  A fixed ceiling
 * compiled into the source cannot do that: the same binary runs on a
 * laptop and on a 500 GB node.
 *
 * "Available" is the smallest of several limits, because any one of them
 * can bind first:
 *
 *   - physical RAM (sysconf, or sysctl on the BSDs/macOS);
 *   - what the kernel thinks is reclaimable right now -- MemAvailable on
 *     Linux, which is what actually decides whether a big allocation
 *     succeeds when other jobs share the node;
 *   - the cgroup memory limit, v2 then v1.  This is the one that binds
 *     under SLURM, Kubernetes, or Docker, where physical RAM is
 *     irrelevant and exceeding the cgroup gets the process OOM-killed
 *     rather than a failed malloc.
 *
 * OQP_MEMORY_LIMIT_GB overrides the probe entirely.  Use it when the
 * automatic answer is wrong -- a batch system that does not use cgroups,
 * or a deliberately smaller budget than the machine allows.
 *
 * Returns 0 when nothing could be determined; callers must treat 0 as
 * "unknown" and skip the check rather than refusing every job.
 */

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if defined(__APPLE__) || defined(__FreeBSD__)
#include <sys/sysctl.h>
#include <sys/types.h>
#endif

#if !defined(_WIN32)
#include <unistd.h>
#endif

/* Smaller of two limits, treating 0 as "no information". */
static uint64_t tighter(uint64_t a, uint64_t b)
{
    if (a == 0) return b;
    if (b == 0) return a;
    return a < b ? a : b;
}

#if defined(__linux__)
/* First value of `key` in /proc/meminfo, in bytes (the file reports kB). */
static uint64_t meminfo_bytes(const char *key)
{
    FILE *f = fopen("/proc/meminfo", "r");
    char line[256];
    size_t keylen = strlen(key);
    uint64_t kb = 0;

    if (!f) return 0;
    while (fgets(line, sizeof line, f)) {
        if (strncmp(line, key, keylen) == 0 && line[keylen] == ':') {
            if (sscanf(line + keylen + 1, "%llu", (unsigned long long *) &kb) != 1)
                kb = 0;
            break;
        }
    }
    fclose(f);
    return kb * 1024ull;
}

/* A cgroup limit file holds either a byte count or "max"/a huge sentinel. */
static uint64_t cgroup_limit(const char *path)
{
    FILE *f = fopen(path, "r");
    char buf[64];
    unsigned long long v;

    if (!f) return 0;
    if (!fgets(buf, sizeof buf, f)) { fclose(f); return 0; }
    fclose(f);

    if (strncmp(buf, "max", 3) == 0) return 0;
    if (sscanf(buf, "%llu", &v) != 1) return 0;
    /* v1 writes a sentinel near UINT64_MAX to mean "unlimited". */
    if (v >= (1ull << 62)) return 0;
    return (uint64_t) v;
}
#endif

uint64_t oqp_available_memory_bytes(void)
{
    uint64_t limit = 0;
    const char *env = getenv("OQP_MEMORY_LIMIT_GB");

    if (env && *env) {
        char *end = NULL;
        double gb = strtod(env, &end);
        if (end != env && gb > 0.0)
            return (uint64_t) (gb * 1073741824.0);
        /* Unparseable: fall through to the probe rather than guessing. */
    }

#if defined(__linux__)
    {
        long pages = sysconf(_SC_PHYS_PAGES);
        long pagesz = sysconf(_SC_PAGESIZE);
        if (pages > 0 && pagesz > 0)
            limit = (uint64_t) pages * (uint64_t) pagesz;
    }
    limit = tighter(limit, meminfo_bytes("MemAvailable"));
    limit = tighter(limit, cgroup_limit("/sys/fs/cgroup/memory.max"));
    limit = tighter(limit,
                    cgroup_limit("/sys/fs/cgroup/memory/memory.limit_in_bytes"));

#elif defined(__APPLE__) || defined(__FreeBSD__)
    {
        uint64_t bytes = 0;
        size_t len = sizeof bytes;
#if defined(HW_MEMSIZE)
        int mib[2] = { CTL_HW, HW_MEMSIZE };
#else
        int mib[2] = { CTL_HW, HW_PHYSMEM };
#endif
        if (sysctl(mib, 2, &bytes, &len, NULL, 0) == 0)
            limit = bytes;
    }

#elif !defined(_WIN32)
    {
        long pages = sysconf(_SC_PHYS_PAGES);
        long pagesz = sysconf(_SC_PAGESIZE);
        if (pages > 0 && pagesz > 0)
            limit = (uint64_t) pages * (uint64_t) pagesz;
    }
#endif

    return limit;
}
