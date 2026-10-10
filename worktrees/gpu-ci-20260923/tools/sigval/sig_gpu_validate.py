#!/usr/bin/env python3
"""MRSF sigma-session validation: GPU engine (src/sigma.cu in libopenqp_gpu.so)
vs the June-validated CPU stub (routec_sig_cpustub.cpp), same v3 ABI, identical
random sessions.  The stub is the referee: it was G-sigma3-gated against native
OpenQP MRSF Davidson to 1.7e-11 Ha in June 2026.

Builds a random symmetric dense B file, random orthogonal Ca/Cb, random
symmetric MO Focks, random trial vectors; runs both libraries for singlet and
triplet kinds and two exchange scales; gates on relative max deviation.

Usage:
  python3 sig_gpu_validate.py --cpulib ./libroutec_sig_cpustub.so \
      --gpulib ../../build/libopenqp_gpu.so [--nbf 30] [--naux 120]
"""
import argparse, ctypes, os, struct, sys, time
import numpy as np


def write_B(path, naux, nbf, rng):
    B = rng.standard_normal((naux, nbf, nbf)) / np.sqrt(naux * nbf)
    B = 0.5 * (B + B.transpose(0, 2, 1))          # each aux slice symmetric
    with open(path, "wb") as f:
        f.write(struct.pack("ii", naux, nbf))     # dense header
        B.astype(np.float64).tofile(f)
    return B


class SigLib:
    def __init__(self, path):
        self.lib = ctypes.CDLL(path, mode=ctypes.RTLD_LOCAL)
        self.lib.routec_sig_init.restype = ctypes.c_int

    def run(self, nbf, Ca, Cb, Fa, Fb, nocca, noccb, kind, scale, bvec, nv, ntrial):
        dp = ctypes.POINTER(ctypes.c_double)
        arr = lambda a: np.asfortranarray(a).ctypes.data_as(dp)
        rc = self.lib.routec_sig_init(
            ctypes.byref(ctypes.c_int(nbf)), arr(Ca), arr(Cb), arr(Fa), arr(Fb),
            ctypes.byref(ctypes.c_int(nocca)), ctypes.byref(ctypes.c_int(noccb)),
            ctypes.byref(ctypes.c_int(kind)))
        if rc != 0:
            raise RuntimeError(f"routec_sig_init rc={rc}")
        self.lib.routec_sig_set_scale(ctypes.byref(ctypes.c_double(scale)))
        sig = np.zeros(ntrial * nv, dtype=np.float64)
        info = ctypes.c_int(-1)
        t0 = time.time()
        self.lib.routec_sig_iter(
            bvec.ctypes.data_as(dp), ctypes.byref(ctypes.c_int(nv)),
            sig.ctypes.data_as(dp), ctypes.byref(info))
        dt = time.time() - t0
        self.lib.routec_sig_free()
        if info.value != 0:
            raise RuntimeError(f"routec_sig_iter info={info.value}")
        return sig, dt


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cpulib", required=True)
    ap.add_argument("--gpulib", required=True)
    ap.add_argument("--nbf", type=int, default=30)
    ap.add_argument("--naux", type=int, default=120)
    ap.add_argument("--nocca", type=int, default=6)
    ap.add_argument("--noccb", type=int, default=4)
    ap.add_argument("--nv", type=int, default=3)
    ap.add_argument("--gate", type=float, default=1e-11)
    ap.add_argument("--seed", type=int, default=7)
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)
    nbf, naux = args.nbf, args.naux
    nocca, noccb, nv = args.nocca, args.noccb, args.nv
    ntrial = nocca * (nbf - noccb)

    bpath = os.path.abspath(f"sigval_B_{naux}x{nbf}.bin")
    write_B(bpath, naux, nbf, rng)
    os.environ["OQP_ROUTEC_B"] = bpath

    Ca = np.linalg.qr(rng.standard_normal((nbf, nbf)))[0]
    Cb = np.linalg.qr(rng.standard_normal((nbf, nbf)))[0]
    Fa = rng.standard_normal((nbf, nbf)); Fa = 0.5 * (Fa + Fa.T)
    Fb = rng.standard_normal((nbf, nbf)); Fb = 0.5 * (Fb + Fb.T)
    bvec = np.asfortranarray(rng.standard_normal((ntrial, nv)))

    cpu = SigLib(args.cpulib)
    gpu = SigLib(args.gpulib)

    print(f"[sigval] nbf={nbf} naux={naux} nocca={nocca} noccb={noccb} "
          f"ntrial={ntrial} nv={nv}")
    worst = 0.0
    for kind in (1, 3):
        for scale in (0.5, 1.0):
            sc, tc = cpu.run(nbf, Ca, Cb, Fa, Fb, nocca, noccb, kind, scale, bvec, nv, ntrial)
            sg, tg = gpu.run(nbf, Ca, Cb, Fa, Fb, nocca, noccb, kind, scale, bvec, nv, ntrial)
            ref = np.max(np.abs(sc))
            rel = np.max(np.abs(sg - sc)) / ref
            worst = max(worst, rel)
            print(f"[sigval] kind={kind} scale={scale}: rel max dev {rel:.3e} "
                  f"(|sigma|max {ref:.3e})  cpu {tc*1e3:.1f} ms  gpu {tg*1e3:.1f} ms")
    verdict = "PASS" if worst <= args.gate else "FAIL"
    print(f"[gate] sigma-session GPU vs CPU-stub: worst rel {worst:.3e}  "
          f"{verdict} (<= {args.gate:g})")
    sys.exit(0 if verdict == "PASS" else 1)


if __name__ == "__main__":
    main()
