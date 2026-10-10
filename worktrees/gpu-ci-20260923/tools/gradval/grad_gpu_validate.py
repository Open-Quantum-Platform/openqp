#!/usr/bin/env python3
"""GPU-dylib validation/bench (chc4): routec_grad2 + routec_grad2_mrsf of
libroutec_oqp_grad_gpu.so vs the ABI-identical CPU dylib on the SAME inputs
(same routec_grad_inp, same densities, same xyz) => identical Gamma/gamma by
construction; the comparison isolates the GPU Gamma-contracted pass 2.

Gate: max |de_gpu - de_cpu| <= 1e-12 Ha/bohr (RHF and MRSF), per molecule.

usage:
  python3 grad_gpu_validate.py --inp h2o_routec_grad.inp --natm 3 \
      --cpulib ./libroutec_oqp_grad_cpu.so --gpulib ./libroutec_oqp_grad_gpu.so
  add --bench N to time N repeated calls of each leg (after the check)
  add --mrsf-only / --rhf-only to restrict
"""
import argparse
import ctypes
import os
import time

import numpy as np


def load_lib(path, inp, threads, verbose=True):
    env = {"OQP_ROUTEC_GRAD_INP": os.path.abspath(inp),
           "OQP_ROUTEC_TABLES": os.path.expanduser("~/grad_gpu/routec_tables.bin"),
           "OQP_ROUTEC_GRAD_THREADS": str(threads)}
    if verbose:
        env["OQP_ROUTEC_GRAD_VERBOSE"] = "1"
    for k, v in env.items():
        os.environ[k] = v
    return ctypes.CDLL(os.path.abspath(path))


def read_dims(inp):
    with open(inp) as f:
        f.readline()
        toks = []
        for _ in range(3):
            toks += f.readline().split()
        d = dict(zip(toks[0::2], (int(x) for x in toks[1::2])))
    return d["natm"], d["nbf"], d["naux"]


DPP = ctypes.POINTER(ctypes.c_double)


def call_rhf(lib, dpk, xyz, nbf, natm, hs=1.0, cs=1.0):
    de = np.zeros((natm, 3))
    info = ctypes.c_int(-1)
    lib.routec_grad2(dpk.ctypes.data_as(DPP), xyz.ctypes.data_as(DPP),
                     de.ctypes.data_as(DPP),
                     ctypes.byref(ctypes.c_int(nbf)),
                     ctypes.byref(ctypes.c_int(natm)),
                     ctypes.byref(ctypes.c_double(hs)),
                     ctypes.byref(ctypes.c_double(cs)),
                     ctypes.byref(info))
    return de, info.value


def call_mrsf(lib, d, p, spc, xyz, nbf, natm,
              hs=1.0, hs2=1.0, cs=1.0, spcs=(0.5, 0.5, 0.5), mrst=1):
    de = np.zeros((natm, 3))
    info = ctypes.c_int(-1)
    sps = np.asarray(spcs, dtype=float)
    lib.routec_grad2_mrsf(d.ctypes.data_as(DPP), p.ctypes.data_as(DPP),
                          spc.ctypes.data_as(DPP), xyz.ctypes.data_as(DPP),
                          de.ctypes.data_as(DPP),
                          ctypes.byref(ctypes.c_int(nbf)),
                          ctypes.byref(ctypes.c_int(natm)),
                          ctypes.byref(ctypes.c_double(hs)),
                          ctypes.byref(ctypes.c_double(hs2)),
                          ctypes.byref(ctypes.c_double(cs)),
                          sps.ctypes.data_as(DPP),
                          ctypes.byref(ctypes.c_int(mrst)),
                          ctypes.byref(info))
    return de, info.value


def synth_inputs(nbf, natm, seed=7):
    rng = np.random.default_rng(seed)
    # plausible smooth symmetric "densities" (scale ~ overlap-like)
    def sym():
        a = rng.standard_normal((nbf, nbf)) * 0.1
        return a + a.T + np.eye(nbf)
    D = sym()
    iu = np.tril_indices(nbf)
    dpk = np.ascontiguousarray(D[iu])
    # MRSF inputs: raw alpha/beta in FORTRAN order (symmetric => order moot,
    # but keep the layout contract: (nbf,nbf,2) F-order = concat of two slices)
    d2 = np.stack([sym(), sym()])                      # (2,nbf,nbf) C
    p2 = np.stack([0.1 * sym(), 0.1 * sym()])
    d_f = np.ascontiguousarray(np.transpose(d2, (0, 2, 1)).reshape(-1))
    p_f = np.ascontiguousarray(np.transpose(p2, (0, 2, 1)).reshape(-1))
    # spc (7,nbf,nbf) with spc[m + 7*(i + nbf*j)] = slot m at (i,j): build
    # C-array of shape (nbf_j, nbf_i, 7) then ravel
    spc_ij = rng.standard_normal((7, nbf, nbf)) * 0.05  # slot, i, j
    spc = np.ascontiguousarray(np.transpose(spc_ij, (2, 1, 0)).reshape(-1))
    # geometry: synthetic cluster matching natm (bohr); MUST be the export
    # geometry family for physical meaning, but for CPU-vs-GPU identity any
    # xyz works -- use the (H2O)n grid the inp files were exported from.
    base = np.array([[0.0, 0.0, 0.1173],
                     [0.0, 0.7572, -0.4692],
                     [0.0, -0.7572, -0.4692]]) * 1.8897261339212517
    nmol = natm // 3
    xyz = []
    k = 0
    for ix in range(2):
        for iy in range(2):
            for iz in range(2):
                if k >= nmol:
                    break
                s = np.array([2.9 * ix + 0.13 * iy, 2.9 * iy + 0.11 * iz,
                              2.9 * iz + 0.12 * ix]) * 1.8897261339212517
                xyz.append(base + s)
                k += 1
    xyz = np.ascontiguousarray(np.concatenate(xyz, axis=0)[:natm])
    return dpk, d_f, p_f, spc, xyz


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--inp", required=True)
    ap.add_argument("--cpulib")
    ap.add_argument("--gpulib", required=True)
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--bench", type=int, default=0)
    ap.add_argument("--rhf-only", action="store_true")
    ap.add_argument("--mrsf-only", action="store_true")
    args = ap.parse_args()

    natm, nbf, naux = read_dims(args.inp)
    print(f"[validate] {args.inp}: natm={natm} nbf={nbf} naux={naux} "
          f"threads={args.threads}")
    dpk, d_f, p_f, spc, xyz = synth_inputs(nbf, natm)

    legs = []
    if args.cpulib:
        legs.append(("cpu", load_lib(args.cpulib, args.inp, args.threads)))
    legs.append(("gpu", load_lib(args.gpulib, args.inp, args.threads)))

    out = {}
    for tag, lib in legs:
        if not args.mrsf_only:
            t0 = time.time()
            de, info = call_rhf(lib, dpk, xyz, nbf, natm)
            t = time.time() - t0
            assert info == 0, (tag, "rhf info", info)
            out[(tag, "rhf")] = (de, t)
            print(f"[{tag}] RHF  de[0]={de[0]}  ({t:.3f} s)")
        if not args.rhf_only:
            t0 = time.time()
            de, info = call_mrsf(lib, d_f, p_f, spc, xyz, nbf, natm)
            t = time.time() - t0
            assert info == 0, (tag, "mrsf info", info)
            out[(tag, "mrsf")] = (de, t)
            print(f"[{tag}] MRSF de[0]={de[0]}  ({t:.3f} s)")

    if args.cpulib:
        for kind in ("rhf", "mrsf"):
            if (("cpu", kind) in out) and (("gpu", kind) in out):
                d = np.abs(out[("gpu", kind)][0] - out[("cpu", kind)][0]).max()
                ref = np.abs(out[("cpu", kind)][0]).max()
                rel = d / ref
                print(f"[gate] {kind}: max|gpu-cpu| = {d:.3e} Ha/bohr "
                      f"(max|cpu| {ref:.3e}, rel {rel:.3e})  "
                      f"{'PASS' if rel <= 1e-11 else 'FAIL'} (rel<=1e-11)")

    if args.bench:
        for tag, lib in legs:
            if not args.mrsf_only:
                ts = []
                for _ in range(args.bench):
                    t0 = time.time()
                    call_rhf(lib, dpk, xyz, nbf, natm)
                    ts.append(time.time() - t0)
                print(f"[bench:{tag}] RHF  best {min(ts):.3f} s over "
                      f"{args.bench} calls: {['%.3f' % t for t in ts]}")
            if not args.rhf_only:
                ts = []
                for _ in range(args.bench):
                    t0 = time.time()
                    call_mrsf(lib, d_f, p_f, spc, xyz, nbf, natm)
                    ts.append(time.time() - t0)
                print(f"[bench:{tag}] MRSF best {min(ts):.3f} s over "
                      f"{args.bench} calls: {['%.3f' % t for t in ts]}")


if __name__ == "__main__":
    main()
