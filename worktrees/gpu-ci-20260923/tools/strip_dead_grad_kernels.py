#!/usr/bin/env python3
"""Strip the DEAD __global__ routec_grad3c_k_* kernels from the generated
derivative .cuh.

The generator emits, per class, three things: grad3c_al_XX_Y and grad3c_lab_XX_Y
(the __forceinline__ __device__ math) AND a standalone __global__
routec_grad3c_k_XX_Y.  grad.cu never launches routec_grad3c_k_* -- it launches
its own gk_XX_Y (DEF_GKERN) which re-inlines grad3c_al/lab.  So every dead
__global__ forces a SECOND full SASS expansion of the al+lab math (both branches
inlined), doubling the device compile with zero runtime benefit.

Removing them roughly halves grad.cu's device codegen and is what makes -O3
tractable on this TU.  Safe: they are unreferenced anywhere in the tree.

Usage: strip_dead_grad_kernels.py <in.cuh> <out.cuh>
"""
import sys, re

src, dst = sys.argv[1], sys.argv[2]
lines = open(src).read().splitlines(keepends=True)
out = []
i = 0
removed = 0
DEAD = re.compile(r'^__global__ void routec_grad3c_k_')
while i < len(lines):
    if DEAD.match(lines[i]):
        # skip from here until braces balance back to zero (top-level '}')
        depth = 0
        started = False
        j = i
        while j < len(lines):
            depth += lines[j].count('{') - lines[j].count('}')
            if '{' in lines[j]:
                started = True
            if started and depth <= 0:
                j += 1
                break
            j += 1
        removed += 1
        i = j
        continue
    out.append(lines[i])
    i += 1

# Also drop the trailing dead dispatch block: the routec_grad3c_cuda_fn typedef,
# the RoutecGradCudaClass struct, the routec_grad3c_cuda[] table (which named the
# kernels above), and nclasses.  grad.cu references none of these.
text = ''.join(out)
m = re.search(r'\ntypedef void \(\*routec_grad3c_cuda_fn\)', text)
if m:
    end = re.search(r'static const int routec_grad3c_cuda_nclasses = \d+;\n', text)
    if end:
        text = text[:m.start()] + '\n' + text[end.end():]
        sys.stderr.write("dropped trailing dead routec_grad3c_cuda[] dispatch table\n")

open(dst, 'w').write(text)
sys.stderr.write(f"stripped {removed} dead __global__ routec_grad3c_k_* kernels; "
                 f"{len(lines)} -> {len(text.splitlines())} lines\n")
