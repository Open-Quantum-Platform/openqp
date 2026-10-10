"""ctypes wrapper for the openqp-gpu library (libopenqp_gpu.so).

Loads the shared library, packs matrices into OpenQP lower-triangular order,
and calls the C SCF drivers.  Point OPENQP_GPU_LIB at the built .so and
OPENQP_GPU_B at the density-fitting tensor.
"""
import os, ctypes, numpy as np

_LIB = None


def _lib():
    global _LIB
    if _LIB is None:
        path = os.environ.get("OPENQP_GPU_LIB", "libopenqp_gpu.so")
        _LIB = ctypes.CDLL(path)
    return _LIB


def _pack(M):
    """Full symmetric (nbf,nbf) -> packed lower-tri t = i*(i+1)/2 + j."""
    nbf = M.shape[0]
    out = np.zeros(nbf * (nbf + 1) // 2)
    for i in range(nbf):
        b = i * (i + 1) // 2
        out[b:b + i + 1] = M[i, :i + 1]
    return np.ascontiguousarray(out)


def solve_rhf(H, S, nocc, enuc, guess="gwh", conv_e=1e-8, conv_d=1e-7, maxit=100,
              scale_coul=1.0, scale_exch=1.0):
    """Closed-shell density-fitting HF/hybrid SCF. Returns (E, cycles).

    H, S: full symmetric (nbf,nbf).  guess: 'gwh' | 'core' | path to packed
    density file.  Set OPENQP_GPU_B to the DF tensor before calling.
    """
    lib = _lib()
    nbf = H.shape[0]
    hpk, spk = _pack(H), _pack(S)
    if guess == "core":
        os.environ["OQP_SCF_COREGUESS"] = "1"
    elif guess != "gwh":
        os.environ["OQP_SCF_GUESS_D"] = guess     # packed density file
    ntri = nbf * (nbf + 1) // 2
    d = np.zeros(ntri); c = np.zeros(nbf * nbf); eps = np.zeros(nbf)
    P = ctypes.POINTER(ctypes.c_double); Pi = ctypes.POINTER(ctypes.c_int)
    dp = lambda x: x.ctypes.data_as(P)
    e = ctypes.c_double(0.0); ncyc = ctypes.c_int(0); info = ctypes.c_int(1)
    ci = lambda v: ctypes.byref(ctypes.c_int(v))
    cd = lambda v: ctypes.byref(ctypes.c_double(v))
    lib.routec_scf_solve(dp(hpk), dp(spk), ci(nbf), ci(nocc), cd(enuc),
                         cd(scale_coul), cd(scale_exch), cd(conv_e), cd(conv_d),
                         ci(maxit), ctypes.byref(e), dp(d), dp(c), dp(eps),
                         ctypes.byref(ncyc), ctypes.byref(info))
    if info.value != 0:
        raise RuntimeError(f"routec_scf_solve failed, info={info.value}")
    return e.value, ncyc.value


def solve_uhf(*a, **k):
    raise NotImplementedError("UHF wrapper: TODO (C entry routec_scf_solve_uhf exists)")
