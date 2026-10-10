#!/usr/bin/env python3
"""COMPLETELY pyscf-free prep for the standalone GPU harness (xcgrid format).

Pipeline (geometry -> full ART dir):
  1. OpenQP MINIMAL invocation (set_basis + ints_1e + guess; NO SCF, ~0.2 s):
     S, Hcore, Hückel guess DM_A, basis (all in OpenQP's CARTESIAN AO frame).
  2. ENUC by direct formula; NOCC from electron count.
  3. Native DF tensor: openqp_gpu_build_df (3c integrals + 2c metric ON GPU,
     ROUTEC_NATIVE_METRIC), perm/dscale transform -> OpenQP frame.
  4. XC grid: self-contained Becke/Treutler/Lebedev generator (grid_native.py,
     Lebedev constants ported from OpenQP's lebedev.F90).
  5. basis.bin / cart.bin in the OpenQP cartesian component order with the
     ANALYTIC normalization bfnrm = sqrt((2l-1)!!/((2a-1)!!(2b-1)!!(2c-1)!!)).
     No c2s: the whole calculation lives in the cartesian frame (ncart == nbf).

Usage: nativeprep_xcgrid.py NW <outdir> [FUNC=BLYP] [TABLES] [BUILDER]
"""
import sys, os, math, struct, subprocess, time, types
# Cap OpenMP threads BEFORE importing oqp / launching the builder. On a big
# many-core box (e.g. 152 cores) an unset OMP_NUM_THREADS defaults to ALL cores;
# under machine load that oversubscribes and the CPU pair/Schwarz build thrashes
# (measured (H2O)16: 7.6s @152thr vs 0.27s @32thr -> 14s vs 4.2s wall). 32 is the
# sweet spot; respect an explicit user setting.
os.environ.setdefault("OMP_NUM_THREADS", str(min(32, os.cpu_count() or 8)))
import numpy as np
import basis_set_exchange as bse
import oqp
# pyoqp's Runner import chain pulls the OPTIONAL geometry-optimizer dep
# libdlfind (runfunc -> libdlfind), which we never use (energy prep only).
# Stub it out if absent so the import succeeds on venvs without dl-find.
try:
    import libdlfind  # noqa: F401
except ImportError:
    def _passthru(*a, **k):
        # works as a plain decorator (@wrapper -> returns the function) and as
        # a decorator factory / arbitrary callable
        if len(a) == 1 and callable(a[0]) and not k: return a[0]
        return _passthru
    class _Stub(types.ModuleType):
        def __getattr__(self, name): return _passthru
    _f = _Stub("libdlfind"); _f.__path__ = []           # package-like
    _cb = _Stub("libdlfind.callback")
    sys.modules["libdlfind"] = _f
    sys.modules["libdlfind.callback"] = _cb
from oqp.pyoqp import Runner
from grid_native import build_grid, write_grid_bin, water_cluster, molecule

t0 = time.time()
NW = int(sys.argv[1]); wd = sys.argv[2]
FUNC = sys.argv[3] if len(sys.argv) > 3 else "BLYP"
TABLES = sys.argv[4] if len(sys.argv) > 4 else "/bighome/cheolho.choi/openqp-gpu/data/routec_tables.bin"
BUILDER = sys.argv[5] if len(sys.argv) > 5 else "/bighome/cheolho.choi/openqp-gpu/build/openqp_gpu_build_df"
HFSCALE = {"BLYP": 0.0, "SLATER": 0.0, "SVWN": 0.0, "B3LYP": 0.2, "BHHLYP": 0.5}[FUNC.upper()]
os.makedirs(wd, exist_ok=True)
# DEFAULT guess = native minao/SAD (validated: bit-identical E, fewer SCF cycles,
# cost hidden behind the B-build). OQP_GUESS=huckel|hcore|sap forces the OpenQP
# guess instead. minao falls back to Hueckel automatically if it can't be built.
GUESS = os.environ.get('OQP_GUESS', 'minao')
# OQPGPU_API=1: call the GPU engine as an in-process LIBRARY (ctypes on
# libopenqp_gpu_df.so) -- shells/H/S pass in MEMORY; no builder subprocess, no
# shells.txt/H.pk/operand files. guess + grid keep their file+poll channels
# (unchanged overlap). Only meaningful with RUN_INPROC.
API = os.environ.get('OQPGPU_API', '1') == '1' and bool(os.environ.get('RUN_INPROC'))

# launch the XC-grid generator IMMEDIATELY (needs only NW) so it overlaps the
# whole OpenQP prep + shells + B build; grid.bin appears atomically (tmp+rename).
GRID_PY = os.environ.get("GRID_PY", sys.executable)
gcmd = [GRID_PY, os.path.join(os.path.dirname(os.path.abspath(__file__)),
        "grid_native.py"), str(NW), f"{wd}/grid.bin"]
pg = subprocess.Popen(gcmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)

# ---- cartesian conventions (same machinery as native_prep.py) ----------------
def _df2(n):
    r = 1
    while n > 1: r *= n; n -= 2
    return r
def _routec_cart(l):
    return [(ax, ay, l-ax-ay) for ax in range(l, -1, -1) for ay in range(l-ax, -1, -1)]
_OQP_CART = {
    0: [(0,0,0)],
    1: [(1,0,0),(0,1,0),(0,0,1)],
    2: [(2,0,0),(0,2,0),(0,0,2),(1,1,0),(1,0,1),(0,1,1)],
    3: [(3,0,0),(0,3,0),(0,0,3),(2,1,0),(2,0,1),(1,2,0),(0,2,1),(1,0,2),(0,1,2),(1,1,1)],
    4: [(4,0,0),(0,4,0),(0,0,4),(3,1,0),(3,0,1),(1,3,0),(0,3,1),(1,0,3),(0,1,3),
        (2,2,0),(2,0,2),(0,2,2),(2,1,1),(1,2,1),(1,1,2)],
}
def cart_transform(l):
    rc = _routec_cart(l); df2l = _df2(2*l-1)
    perm = [rc.index(t) for t in _OQP_CART[l]]
    dsc  = [math.sqrt(df2l / (_df2(2*a-1)*_df2(2*b-1)*_df2(2*c-1)))
            for (a,b,c) in _OQP_CART[l]]
    return perm, dsc

# ---- 1. OpenQP minimal invocation (NO SCF) ------------------------------------
atoms = molecule(NW)                            # (H2O)NW, or OQP_SYSTEM_XYZ molecule
ZNUM = {"H":1,"HE":2,"LI":3,"BE":4,"B":5,"C":6,"N":7,"O":8,"F":9,"NE":10,
        "NA":11,"MG":12,"AL":13,"SI":14,"P":15,"S":16,"CL":17,"AR":18}
SYSTEM = "".join(f"\n   {ZNUM[s.upper()]}   {x:.9f} {y:.9f} {z:.9f}" for s, (x, y, z) in atoms)
CFG = {'input': {'system': SYSTEM, 'charge': '0', 'runtype': 'energy',
                 'basis': 'cc-pvdz', 'method': 'hf', 'functional': '', 'd4': 'False'},
       # 'minao' is OUR native SAD guess (native_minao.py, applied after the
       # shell block); OpenQP runs a throwaway Hueckel we then overwrite.
       'guess': {'type': ('huckel' if GUESS == 'minao' else GUESS)},
       'scf':   {'type': 'rhf', 'multiplicity': '1',
                 'save_molden': 'False', 'incremental': 'False'}}
if os.environ.get('OQP_ISPHER'):                 # newer OpenQP defaults cc-pVDZ to
    CFG['input']['ispher'] = os.environ['OQP_ISPHER']  # spherical; force cartesian
r = Runner(project='nprep', input_dict=CFG, log=f'{wd}/nprep.log', silent=1, usempi=False)
mol = r.mol
try:                                            # build-version guard (see agent report)
    mol.data['OQP::log_filename']
except AttributeError:
    mol.data['OQP::log_filename'] = mol.log
    oqp.oqp_banner(mol)
oqp.library.set_basis(mol)
oqp.library.ints_1e(mol)                        # -> OQP::SM, OQP::TM, OQP::Hcore
oqp.library.guess(mol)                          # -> OQP::DM_A (Hückel)
t_oqp = time.time()

Spk = np.array(mol.data['OQP::SM'])
Hpk = np.array(mol.data['OQP::Hcore'])
DA  = np.array(mol.data['OQP::DM_A'])
bas = mol.data.get_basis()
xyz = np.asarray(mol.get_system(), float).reshape(-1, 3)     # Bohr
Zs  = np.array([ZNUM[s.upper()] for s, _ in atoms], float)
nbf = int((math.isqrt(8*len(Spk)+1)-1)//2)
ntri = nbf*(nbf+1)//2
assert ntri == len(Spk)

def unpack(p):
    M = np.zeros((nbf, nbf))
    iu, ju = np.tril_indices(nbf)
    M[iu, ju] = p; M[ju, iu] = p
    return M
S = unpack(Spk); H = unpack(Hpk)
if not API:   # library mode passes H/S in memory
    np.save(f"{wd}/H_w{NW}.npy", H); np.save(f"{wd}/S_w{NW}.npy", S)
    Hpk.astype("<f8").tofile(f"{wd}/H.pk"); Spk.astype("<f8").tofile(f"{wd}/S.pk")  # raw packed for the in-process solver

# guess: the solver expects the TOTAL (occ-2) density; OpenQP RHF DM_A may be
# alpha-only. Detect via Tr(D S) and rescale.
nelec = int(Zs.sum())
# packed-tri trace: Tr(D S) = sum_diag D_ii S_ii + 2 sum_offdiag D_ij S_ij
diag_idx = np.array([i*(i+1)//2 + i for i in range(nbf)])
w2 = np.full(ntri, 2.0); w2[diag_idx] = 1.0
trDS = float(np.sum(DA * Spk * w2))
scale = nelec / trDS
# minao mode writes its own guess later (overlapped with the builder); the
# builder polls for the file, so DON'T pre-write the Hueckel one here or the
# poll would race and read Hueckel.
if GUESS != 'minao':
    (DA * scale).astype("<f8").tofile(f"{wd}/guess_w{NW}.bin")

iu2, ju2 = np.triu_indices(len(Zs), 1)
enuc = float(np.sum(Zs[iu2]*Zs[ju2] / np.linalg.norm(xyz[iu2]-xyz[ju2], axis=1)))

# ---- 2. shells + operands for the GPU builder (native_prep machinery) --------
angs = np.asarray(bas['angs']); ncontr = np.asarray(bas['ncontr'])
centers = np.asarray(bas['centers']); alpha = np.asarray(bas['alpha']); coef = np.asarray(bas['coef'])
lines = []; nsh = 0; ip = 0; nao = 0
sh_am = []; sh_g0 = []; sh_nc = []; sh_ao = []; sh_cx = []; sh_cy = []; sh_cz = []
ex_all = []; cc_all = []
for ish in range(len(angs)):
    l = int(angs[ish]); nc = int(ncontr[ish])
    e = alpha[ip:ip+nc]; c = coef[ip:ip+nc]; ip += nc
    x, y, z = xyz[int(centers[ish])]
    lines.append(f"{l} {x:.12f} {y:.12f} {z:.12f} {nc}")
    for a, cc in zip(e, c): lines.append(f"{a:.12g} {cc:.12g}")
    sh_am.append(l); sh_g0.append(len(ex_all)+1); sh_nc.append(nc); sh_ao.append(nao+1)
    sh_cx.append(x); sh_cy.append(y); sh_cz.append(z)
    ex_all += list(e); cc_all += list(c)
    nsh += 1; nao += (l+1)*(l+2)//2
if not API:
    open(f"{wd}/shells_c_w{NW}.txt", "w").write(f"{nsh}\n" + "\n".join(lines) + "\n")
assert nao == nbf, f"cartesian count mismatch {nao} vs {nbf} (OpenQP not cartesian?)"

def gto_norm(l, a):
    return math.sqrt(2**(2*l+3)*math.factorial(l+1)*(2*a)**(l+1.5)
                     /(math.factorial(2*l+2)*math.sqrt(math.pi)))
auxcache = {}; alines = []; ansh = 0; naux = 0
def aux_shells_for(znum):
    if znum in auxcache: return auxcache[znum]
    el = list(bse.get_basis("def2-universal-jkfit", elements=[int(znum)])["elements"].values())[0]
    out = []
    for sh in el["electron_shells"]:
        exps = np.array([float(x) for x in sh["exponents"]])
        for l, crow in zip(sh["angular_momentum"]*len(sh["coefficients"]), sh["coefficients"]):
            c = np.array([float(x) for x in crow]); m = np.abs(c) > 0
            e = exps[m]; cn = c[m]*np.array([gto_norm(l, a) for a in e])
            out.append((l, e, cn))
    auxcache[znum] = out; return out
a_am = []; a_x = []; a_y = []; a_z = []; a_nc = []; a_ex = []; a_cc = []
for ia, z in enumerate(Zs):
    x, y, zz = xyz[ia]
    for l, e, cn in aux_shells_for(z):
        alines.append(f"{l} {x:.12f} {y:.12f} {zz:.12f} {len(e)}")
        for a, cc in zip(e, cn): alines.append(f"{a:.12g} {cc:.12g}")
        a_am.append(l); a_x.append(x); a_y.append(y); a_z.append(zz); a_nc.append(len(e))
        a_ex += list(e); a_cc += list(cn)
        ansh += 1; naux += 2*l+1
if not API:
    open(f"{wd}/aux_c_w{NW}.txt", "w").write(f"{ansh}\n" + "\n".join(alines) + "\n")

npair = nao*(nao+1)//2
perm = np.arange(nao, dtype=np.int32); dscale = np.ones(nao); o = 0
for ish in range(nsh):
    l = sh_am[ish]; ncmp = (l+1)*(l+2)//2
    pm, ds = cart_transform(l)
    for b in range(ncmp): perm[o+b] = o + pm[b]; dscale[o+b] = ds[b]
    o += ncmp
if not API:   # library mode derives every transform operand internally
    ilo, jlo = np.tril_indices(nao)
    np.array([nao, naux, npair], dtype=np.int32).tofile(f"{wd}/meta_w{NW}.bin")
    np.zeros((naux, naux)).tofile(f"{wd}/cag_w{NW}.bin")
    np.zeros((naux, naux)).tofile(f"{wd}/linv_w{NW}.bin")
    perm.tofile(f"{wd}/perm_w{NW}.bin")
    np.ascontiguousarray(dscale, dtype=np.float64).tofile(f"{wd}/dscale_w{NW}.bin")
    np.arange(naux, dtype=np.int32).tofile(f"{wd}/aperm_w{NW}.bin")
    np.ones(naux, dtype=np.float64).tofile(f"{wd}/ascale_w{NW}.bin")
    ilo.astype(np.int32).tofile(f"{wd}/iu_w{NW}.bin")
    jlo.astype(np.int32).tofile(f"{wd}/ju_w{NW}.bin")

# ---- native minao/SAD guess (deferred; runs CONCURRENTLY with the builder) ----
# Projected superposition of atomic densities -- a better SCF starting point than
# Hueckel (DFT (H2O)16: 10 SCF cycles vs 16). Its analytic-overlap build is a
# couple of CPU-seconds, so we hide it behind the builder's GPU B-build: launch
# the builder first, compute minao here, write the guess ATOMICALLY, and the
# builder polls for the file just before its (post-B-build) solve.
def write_minao_guess():
    _tm = time.time()
    tmp = f"{wd}/guess_w{NW}.bin.tmp"
    try:
        from native_minao import minao_guess
        D0 = minao_guess(sh_am, sh_cx, sh_cy, sh_cz, sh_g0, sh_nc,
                         ex_all, cc_all, dscale, S, Zs, xyz, verbose=True)
        d0pk = D0[np.tril_indices(nbf)]             # OpenQP row-walk packed order
        trDS0 = float(np.sum(d0pk * Spk * w2))
        (d0pk * (nelec / trDS0)).astype("<f8").tofile(tmp)
        print(f"[minao] compute {time.time()-_tm:.2f}s trDS={trDS0:.4f} rescale->{nelec}", flush=True)
    except Exception as e:                          # never break a run: fall back to Hueckel
        print(f"[minao] FAILED ({type(e).__name__}: {e}); falling back to Hueckel", flush=True)
        (DA * scale).astype("<f8").tofile(tmp)
    os.replace(tmp, f"{wd}/guess_w{NW}.bin")        # atomic: builder never sees a partial file
t_shell = time.time()

# ---- 5. basis.bin / cart.bin (OpenQP cartesian order; analytic bfnrm) ---------
bfnrm = dscale.copy()                    # AO_oqp = dscale * raw collocation
maxang = max(sh_am); maxcart = (maxang+1)*(maxang+2)//2
def wi(f, *v): f.write(struct.pack("<%dq" % len(v), *v))
def wdv(f, a): f.write(np.ascontiguousarray(a, dtype="<f8").tobytes())
with open(f"{wd}/basis.bin", "wb") as f:
    wi(f, nsh, nao, len(ex_all))
    for s in range(nsh):
        wi(f, sh_am[s], 0, sh_g0[s], sh_nc[s], sh_ao[s], 0)
        wdv(f, [sh_cx[s], sh_cy[s], sh_cz[s]]); wdv(f, [1e30])
    wdv(f, ex_all); wdv(f, cc_all); wdv(f, [1e30]*len(ex_all)); wdv(f, bfnrm)
tx = np.zeros(maxcart*(maxang+1), dtype=np.int64); ty = tx.copy(); tz = tx.copy()
for l in range(maxang+1):
    for k, (a, b, c) in enumerate(_OQP_CART[l]):
        tx[k+maxcart*l] = a; ty[k+maxcart*l] = b; tz[k+maxcart*l] = c
with open(f"{wd}/cart.bin", "wb") as f:
    wi(f, maxang, maxcart); f.write(tx.tobytes()); f.write(ty.tobytes()); f.write(tz.tobytes())


# ---- 3+4. native GPU B build CONCURRENT with the XC grid subprocess -----------
# (independent outputs; the grid's small kernels coexist fine with the builder)
env = dict(os.environ); env.update(ROUTEC_CART="1", ROUTEC_NATIVE_METRIC="1", ROUTEC_B_SIGMA="1")
cmd = [BUILDER, TABLES, f"{wd}/shells_c_w{NW}.txt", f"{wd}/aux_c_w{NW}.txt",
       f"{wd}/M_scratch.bin", "1e-13", wd, f"w{NW}", "0", f"{wd}/B_w{NW}.bin"]  # nrep=0: production (no bench)
GRID_PY = os.environ.get("GRID_PY", sys.executable)
gcmd = [GRID_PY, os.path.join(os.path.dirname(os.path.abspath(__file__)),
        "grid_native.py"), str(NW), f"{wd}/grid.bin"]
RUN = os.environ.get("RUN_INPROC", "")          # "HF"/"DFT": solve inside the builder
if RUN:
    cmd = cmd[:9]                               # drop Bout: no B file at all
    env.update(ROUTEC_SOLVE_HS=f"{wd}/", ROUTEC_SOLVE_NOCC=str(nelec//2),
               ROUTEC_SOLVE_ENUC=f"{enuc:.12f}",
               ROUTEC_SOLVE_SE=("1.0" if RUN=="HF" else f"{HFSCALE:.6f}"),
               OQP_SCF_GUESS_D=f"{wd}/guess_w{NW}.bin")
    if RUN=="DFT":
        env.update(OQP_OWNXC_DIR=wd, OQP_OWNXC_FUNC=FUNC)
    # launch the engine FIRST, then compute the minao guess concurrently with its
    # B-build (the engine polls for guess_w*.bin before the solve). For Hueckel
    # the guess is already on disk, so this just runs the engine.
    if API:
        # LIBRARY MODE: one in-process ctypes call; shells/H/S pass in MEMORY.
        # CDLL releases the GIL, so a plain thread overlaps the GPU build+solve
        # with the minao construction exactly like the subprocess did.
        for k, v in env.items():
            if k not in os.environ or os.environ[k] != v: os.environ[k] = v
        import threading, oqpgpu
        libdf = os.path.join(os.path.dirname(BUILDER), "libopenqp_gpu_df.so")
        argvL = [BUILDER, TABLES, "", "", f"{wd}/M_scratch.bin", "1e-13",
                 wd, f"w{NW}", "0"]
        orb  = (sh_am, sh_cx, sh_cy, sh_cz, sh_nc, ex_all, cc_all)
        auxs = (a_am, a_x, a_y, a_z, a_nc, a_ex, a_cc)
        holder = {}
        def _call():
            try: holder["rc"] = oqpgpu.run(libdf, argvL, orb, auxs, Hpk, Spk)
            except Exception as e: holder["rc"] = -1; holder["err"] = repr(e)
        th = threading.Thread(target=_call); th.start()
        if GUESS == 'minao':
            write_minao_guess()
        th.join()
        if holder.get("rc", -1) != 0:
            print(holder.get("err", ""))
            pg.kill()
            raise SystemExit(f"in-process (library) solve failed rc={holder.get('rc')}")
    else:
        proc = subprocess.Popen(cmd, env=env, stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE, text=True)
        if GUESS == 'minao':
            write_minao_guess()
        out, errout = proc.communicate()
        rc = types.SimpleNamespace(returncode=proc.returncode, stdout=out, stderr=errout)
        for ln in rc.stdout.splitlines():
            if "inproc" in ln: print(ln, flush=True)
        if rc.returncode != 0 or "inproc-solve" not in rc.stdout:
            print(rc.stdout[-3000:]); print(rc.stderr[-3000:])
            raise SystemExit(f"in-process solve failed rc={rc.returncode}")
    gout, _ = pg.communicate()
    if pg.returncode != 0:
        print(gout[-2000:]); raise SystemExit("grid generation failed")
else:
    # B-to-disk path (no in-process solve): write the minao guess for a later
    # offline solve; no concurrency benefit here so just do it before the build.
    if GUESS == 'minao':
        write_minao_guess()
    rc = subprocess.run(cmd, env=env, capture_output=True, text=True)
    if rc.returncode != 0 or not os.path.exists(f"{wd}/B_w{NW}.bin"):
        print(rc.stdout[-3000:]); print(rc.stderr[-3000:])
        pg.kill()
        raise SystemExit(f"builder failed rc={rc.returncode}")
t_build = time.time()
if not RUN:
    gout, _ = pg.communicate()
    if pg.returncode != 0:
        print(gout[-2000:]); raise SystemExit("grid generation failed")
W = np.fromfile(f"{wd}/grid.bin", dtype=np.int64, count=1)  # npts for the log
t_grid = time.time()

with open(f"{wd}/meta.env", "w") as f:
    f.write(f"NAO={nao}\nNOCC={nelec//2}\nENUC={enuc:.12f}\nNAUX={naux}\n"
            f"HFSCALE={HFSCALE:.6f}\nEREF=0.0\nFUNC={FUNC}\nNCART={nao}\n")
t_end = time.time()
print(f"[nativeprep w{NW}] nao(cart)={nao} naux={naux} npts={int(W[0])} "
      f"trDS={trDS:.3f}->scale{scale:.3f} | oqp {t_oqp-t0:.1f}s shells {t_shell-t_oqp:.1f}s "
      f"Bbuild {t_build-t_shell:.1f}s grid {t_grid-t_build:.1f}s io {t_end-t_grid:.1f}s "
      f"TOTAL {t_end-t0:.1f}s", flush=True)
