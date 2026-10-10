"""ctypes bridge to libopenqp_gpu_df.so: the GPU engine as a standalone library.

run() hands the orbital/aux shells and the packed H/S to oqpgpu_run() IN MEMORY
-- no builder subprocess, no shells.txt / H.pk / operand files. The small async
artifacts keep their file+poll channels so the caller can overlap their
construction with the GPU B-build: write the guess density to the path in
OQP_SCF_GUESS_D (atomically) and, for DFT, grid.bin under OQP_OWNXC_DIR.
Knobs travel via os.environ (getenv in the lib sees them) and an argv list:

  argv = [name, tables.bin, "", "", scratch, tau, prep_dir, tag, nrep]

(argv[2]/[3] shell files are ignored when memory shells are given.)
CDLL calls release the GIL, so calling run() from a threading.Thread overlaps
the whole GPU build+solve with python-side work (minao guess, grid wait).
"""
import ctypes as ct
import numpy as np

class _Shells(ct.Structure):
    _fields_ = [("nsh", ct.c_int),
                ("l",   ct.POINTER(ct.c_int)),
                ("x",   ct.POINTER(ct.c_double)),
                ("y",   ct.POINTER(ct.c_double)),
                ("z",   ct.POINTER(ct.c_double)),
                ("np",  ct.POINTER(ct.c_int)),
                ("ex",  ct.POINTER(ct.c_double)),
                ("cc",  ct.POINTER(ct.c_double))]

class _MemIn(ct.Structure):
    _fields_ = [("orb", _Shells), ("aux", _Shells),
                ("Hpk", ct.POINTER(ct.c_double)),
                ("Spk", ct.POINTER(ct.c_double))]

def _mk_shells(sh, keep):
    """sh = (am, x, y, z, nprim, ex, cc) arrays/lists -> _Shells (keep owns memory)."""
    am, x, y, z, npr, ex, cc = sh
    am  = np.ascontiguousarray(am,  np.int32)
    x   = np.ascontiguousarray(x,   np.float64)
    y   = np.ascontiguousarray(y,   np.float64)
    z   = np.ascontiguousarray(z,   np.float64)
    npr = np.ascontiguousarray(npr, np.int32)
    ex  = np.ascontiguousarray(ex,  np.float64)
    cc  = np.ascontiguousarray(cc,  np.float64)
    assert ex.size == cc.size == int(npr.sum()), "primitive count mismatch"
    keep += [am, x, y, z, npr, ex, cc]
    s = _Shells()
    s.nsh = len(am)
    s.l   = am.ctypes.data_as(ct.POINTER(ct.c_int))
    s.x   = x.ctypes.data_as(ct.POINTER(ct.c_double))
    s.y   = y.ctypes.data_as(ct.POINTER(ct.c_double))
    s.z   = z.ctypes.data_as(ct.POINTER(ct.c_double))
    s.np  = npr.ctypes.data_as(ct.POINTER(ct.c_int))
    s.ex  = ex.ctypes.data_as(ct.POINTER(ct.c_double))
    s.cc  = cc.ctypes.data_as(ct.POINTER(ct.c_double))
    return s

def run(libpath, argv, orb, aux, Hpk, Spk):
    """One in-process geometry->energy call. Returns the lib's exit code (0=ok).
    orb/aux: (am, x, y, z, nprim, ex, cc); Hpk/Spk: packed lower-tri float64."""
    lib = ct.CDLL(libpath)
    lib.oqpgpu_run.restype = ct.c_int
    lib.oqpgpu_run.argtypes = [ct.c_int, ct.POINTER(ct.c_char_p), ct.POINTER(_MemIn)]
    keep = []
    mem = _MemIn()
    mem.orb = _mk_shells(orb, keep)
    mem.aux = _mk_shells(aux, keep)
    Hpk = np.ascontiguousarray(Hpk, np.float64)
    Spk = np.ascontiguousarray(Spk, np.float64)
    keep += [Hpk, Spk]
    mem.Hpk = Hpk.ctypes.data_as(ct.POINTER(ct.c_double))
    mem.Spk = Spk.ctypes.data_as(ct.POINTER(ct.c_double))
    av = (ct.c_char_p * len(argv))(*[a.encode() for a in argv])
    return lib.oqpgpu_run(len(argv), av, ct.byref(mem))
