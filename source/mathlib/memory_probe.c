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

/* What a cgroup directory still allows: its limit minus what is already
 * charged against it.  A cap is not an allowance -- half of a 16 GB cgroup
 * may already be resident (this process's SCF data, or sibling MPI ranks in
 * the same job cgroup), and a guard that compares a new allocation against
 * the cap admits jobs the kernel will kill.  Clamped to 1 byte when the
 * charge has reached the cap, because 0 means "no information" to tighter().
 * No limit file, or an unlimited one, yields 0: nothing to subtract from. */
static uint64_t cgroup_remaining(const char *dir, const char *limit_leaf,
                                 const char *usage_leaf)
{
    char file[4352];
    uint64_t lim, use;

    if (snprintf(file, sizeof file, "%s/%s", dir, limit_leaf)
            >= (int) sizeof file)
        return 0;
    lim = cgroup_limit(file);
    if (lim == 0) return 0;

    use = 0;
    if (snprintf(file, sizeof file, "%s/%s", dir, usage_leaf)
            < (int) sizeof file)
        use = cgroup_limit(file);   /* plain byte count; same parser */

    return use < lim ? lim - use : 1;
}

/* Tightest remaining allowance on the path from this process's own cgroup up
 * to `mount`.  The effective limit of a nested cgroup is the minimum over its
 * ancestors, and under Slurm or systemd the job's cgroup IS nested -- probing
 * only the mount root reads the (usually absent) top-level limit and misses
 * the allocation that actually binds.  Inside a cgroup namespace, the common
 * container case, /proc/self/cgroup reports "/" and this walk degenerates to
 * exactly the root probe it generalises. */
static uint64_t cgroup_walk(const char *mount, const char *relpath,
                            const char *limit_leaf, const char *usage_leaf)
{
    char path[4096];
    uint64_t best = 0;
    size_t base;

    if (snprintf(path, sizeof path, "%s%s", mount, relpath)
            >= (int) sizeof path)
        return 0;
    base = strlen(mount);

    for (;;) {
        char *cut;
        best = tighter(best, cgroup_remaining(path, limit_leaf, usage_leaf));
        cut = strrchr(path, '/');
        if (!cut || (size_t)(cut - path) < base) break;
        *cut = '\0';
    }
    return best;
}

/* The process's own cgroup limit, resolved through /proc/self/cgroup.
 * v2 lines read "0::<path>"; v1 lists one line per controller and the
 * memory one reads "<n>:memory:<path>". */
static uint64_t cgroup_self_limit(void)
{
    FILE *f = fopen("/proc/self/cgroup", "r");
    char line[4096];
    uint64_t best = 0;

    if (!f) return 0;
    while (fgets(line, sizeof line, f)) {
        char *nl = strchr(line, '\n');
        if (nl) *nl = '\0';
        if (strncmp(line, "0::", 3) == 0) {
            best = tighter(best,
                cgroup_walk("/sys/fs/cgroup", line + 3,
                            "memory.max", "memory.current"));
        } else {
            char *c1 = strchr(line, ':');
            if (!c1) continue;
            char *c2 = strchr(c1 + 1, ':');
            if (!c2) continue;
            *c2 = '\0';
            /* the controller field may list several, comma separated */
            if (strstr(c1 + 1, "memory")) {
                best = tighter(best,
                    cgroup_walk("/sys/fs/cgroup/memory", c2 + 1,
                                "memory.limit_in_bytes",
                                "memory.usage_in_bytes"));
            }
        }
    }
    fclose(f);
    return best;
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
    /* The job's own (possibly nested) cgroup first; the mount-root probes
     * stay as a fallback for setups where /proc/self/cgroup is unreadable. */
    limit = tighter(limit, cgroup_self_limit());
    limit = tighter(limit, cgroup_remaining("/sys/fs/cgroup",
                                            "memory.max", "memory.current"));
    limit = tighter(limit, cgroup_remaining("/sys/fs/cgroup/memory",
                                            "memory.limit_in_bytes",
                                            "memory.usage_in_bytes"));

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
