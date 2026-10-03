"""MRSF spin-orbit coupling: the native H_SOC against an independent determinant oracle.

``soc_mrsf`` (source/modules/soc_mrsf.F90) builds the spin-component transition
densities of the MRSF states (compute_tdm) and contracts them with the MO SOC
integrals (compute_soc_matrix).  This test rebuilds H_SOC without that code:

* every MRSF state is expanded in Slater determinants represented as sorted
  tuples of occupied spin orbitals (textbook phases (-1)**(occupied before
  the target));
* the triplet sublevels are |T,+-1> = S+-|T,0>/sqrt2 (no Wigner-Eckart
  assumption);
* the spin-orbital transition densities D[P,Q] = <bra|a+_P a_Q|ket> are
  contracted with the MO integrals the native run exports (OQP::soc_lmo_1e,
  OQP::soc_lmo_2e; l_b real antisymmetric, L_b = -i l_b);
* the singlet-triplet block uses the native Wigner-Eckart form from <S|D|T,0>
  (D_bb = -D_aa), the triplet-triplet block the explicit sublevels.

The native H_SOC (OQP::soc_hsoc_re/_im, a.u.) must agree to 1e-8 cm-1.  Pure
numpy; no external quantum-chemistry package.  Each case runs in a fresh
interpreter (OpenQP keeps process-global state); skipped when liboqp is absent.
"""
import json
import math
import os
import subprocess
import sys
import tempfile
import textwrap
import unittest
from collections import defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SQ2 = math.sqrt(2.0)
FINE_STRUCTURE = 7.2973525693e-3           # source/physical_constants.F90
HA_TO_WAVENUM = 219474.6313708
DFAC = FINE_STRUCTURE ** 2 / 2.0 * HA_TO_WAVENUM


# ---------------------------------------------------------------------------
# Determinant algebra on sorted occupation tuples
# ---------------------------------------------------------------------------

def annihilate(det, k):
    if k not in det:
        return 0, None
    pos = det.index(k)
    return (-1 if pos % 2 else 1), det[:pos] + det[pos + 1:]


def create(det, k):
    if k in det:
        return 0, None
    pos = sum(1 for q in det if q < k)
    new = list(det)
    new.insert(pos, k)
    return (-1 if pos % 2 else 1), tuple(new)


def apply_ops(state, ops, coef=1.0):
    """ops: [('c'|'a', k), ...] written left to right as in the operator product."""
    out = defaultdict(float)
    for det, c in state.items():
        phase, cur = 1, det
        for kind, k in reversed(ops):
            p, cur = (annihilate if kind == 'a' else create)(cur, k)
            if p == 0:
                break
            phase *= p
        else:
            out[cur] += coef * phase * c
    return out


def add_into(acc, state, coef=1.0):
    for d, c in state.items():
        acc[d] += coef * c
    return acc


def dot(s1, s2):
    if len(s1) > len(s2):
        s1, s2 = s2, s1
    return sum(c * s2.get(d, 0.0) for d, c in s1.items())


class SpinOrbitalModel:
    """Spin orbital p alpha -> p, p beta -> norb + p."""

    def __init__(self, norb):
        self.norb = norb

    def a(self, p):
        return p

    def b(self, p):
        return self.norb + p

    def vacuum_string(self, creations):
        det, phase = (), 1
        for k in reversed(creations):
            p, det = create(det, k)
            phase *= p
        return {det: float(phase)}

    def one_body(self, state, t, kind):
        sp = {'a': self.a, 'b': self.b}
        cs, ans = sp[kind[0]], sp[kind[1]]
        out = defaultdict(float)
        for p, q in np.argwhere(np.abs(t) > 0):
            add_into(out, apply_ops(state, [('c', cs(p)), ('a', ans(q))]), t[p, q])
        return out

    def s_plus(self, state):
        return self.one_body(state, np.eye(self.norb), 'ab')

    def s_minus(self, state):
        return self.one_body(state, np.eye(self.norb), 'ba')

    def tdm(self, bra, ket, kind):
        """D[p,q] = <bra| a+_{p s} a_{q s'} |ket> for kind in aa, bb, ab, ba."""
        n = self.norb
        sp = {'a': self.a, 'b': self.b}
        cs, ans = sp[kind[0]], sp[kind[1]]
        d = np.zeros((n, n))
        # only annihilators present in the ket and creators present in the bra matter
        qs = sorted({k for det in ket for k in det})
        ps = sorted({k for det in bra for k in det})
        for q in range(n):
            if ans(q) not in qs:
                continue
            for p in range(n):
                if cs(p) not in ps:
                    continue
                v = apply_ops(ket, [('c', cs(p)), ('a', ans(q))])
                if v:
                    d[p, q] = dot(bra, v)
        return d


class MRSFMap:
    """Packed MRSF vector -> determinant expansion (two-SOMO, M_S = 0)."""

    def __init__(self, norb, nocca, noccb):
        self.m = SpinOrbitalModel(norb)
        self.norb, self.nocca, self.noccb = norb, nocca, noccb
        self.o1, self.o2 = nocca - 2, nocca - 1
        m = self.m
        core = []
        for c in range(noccb):
            core += [m.a(c), m.b(c)]
        self.rplus = m.vacuum_string(core + [m.a(self.o1), m.a(self.o2)])
        self.rminus = m.vacuum_string(core + [m.b(self.o1), m.b(self.o2)])

    def eplus(self, i, a):
        return apply_ops(self.rplus, [('c', self.m.b(a)), ('a', self.m.a(i))])

    def eminus(self, i, a):
        return apply_ops(self.rminus, [('c', self.m.a(a)), ('a', self.m.b(i))])

    def state(self, x, mult):
        lam = -1.0 if mult == 1 else 1.0
        o1, o2 = self.o1, self.o2
        oo = {(o1, o1), (o2, o2), (o2, o1), (o1, o2)}
        x = np.asarray(x).reshape(self.norb - self.noccb, self.nocca)   # [a - noccb, i]
        psi = defaultdict(float)
        for ia in range(x.shape[0]):
            a = ia + self.noccb
            for i in range(self.nocca):
                c = x[ia, i]
                if c == 0.0 or (i, a) in oo:
                    continue
                add_into(psi, self.eplus(i, a), c / SQ2)
                add_into(psi, self.eminus(i, a), lam * c / SQ2)
        cl = x[o1 - self.noccb, o1]
        add_into(psi, self.eplus(o1, o1), cl / SQ2)
        add_into(psi, self.eplus(o2, o2), lam * cl / SQ2)
        if mult == 1:
            add_into(psi, self.eplus(o2, o1), x[o1 - self.noccb, o2])   # G
            add_into(psi, self.eplus(o1, o2), x[o2 - self.noccb, o1])   # D
        return {d: c for d, c in psi.items() if abs(c) > 1e-15}

    def ladder(self, psi0):
        """|T,+-1> = S+-|T,0>/sqrt2."""
        up, dn = self.m.s_plus(psi0), self.m.s_minus(psi0)
        return {-1: {d: c / SQ2 for d, c in dn.items() if abs(c) > 1e-15},
                0: psi0,
                1: {d: c / SQ2 for d, c in up.items() if abs(c) > 1e-15}}


# ---------------------------------------------------------------------------
# SOC assembly from the oracle densities (native conventions)
# ---------------------------------------------------------------------------

def soc_element(model, bra, ket, lmo):
    """<bra| sum_i L.s |ket> in a.u. (no alpha^2/2); l_b(t,u) pairs with D[u,t]."""
    lx, ly, lz = lmo
    daa, dbb = model.tdm(bra, ket, 'aa').T, model.tdm(bra, ket, 'bb').T
    dab, dba = model.tdm(bra, ket, 'ab').T, model.tdm(bra, ket, 'ba').T
    celm_aa, celm_bb = -0.5j * lz, 0.5j * lz
    celm_ba, celm_ab = 0.5 * (ly - 1j * lx), 0.5 * (-ly - 1j * lx)
    return np.sum(celm_aa * daa + celm_bb * dbb + celm_ba * dba + celm_ab * dab)


def oracle_hsoc(xs, xt, nbf, nocca, noccb, lmo):
    ns, nt = len(xs), len(xt)
    mp = MRSFMap(nbf, nocca, noccb)
    S = [mp.state(x, 1) for x in xs]
    T = [mp.ladder(mp.state(x, 3)) for x in xt]
    lx, ly, lz = lmo
    n = ns + 3 * nt
    h = np.zeros((n, n), dtype=complex)
    for i in range(ns):
        for j in range(nt):
            daa = mp.m.tdm(S[i], T[j][0], 'aa').T
            h[i, ns + 3 * j + 1] = np.sum(2.0 * (-0.5j * lz) * daa)
            h[i, ns + 3 * j + 2] = np.sum(0.5 * (ly - 1j * lx) * (-SQ2) * daa)
            h[i, ns + 3 * j + 0] = np.sum(0.5 * (-ly - 1j * lx) * (+SQ2) * daa)
    for i in range(nt):
        for j in range(nt):
            for a in (-1, 0, 1):
                for b in (-1, 0, 1):
                    h[ns + 3 * i + a + 1, ns + 3 * j + b + 1] = soc_element(mp.m, T[i][a], T[j][b], lmo)
    for i in range(ns):
        for j in range(ns, n):
            h[j, i] = np.conj(h[i, j])
    return h


DRIVER = textwrap.dedent('''
    import json, math, sys
    import numpy as np
    from oqp.pyoqp import Runner
    functional = sys.argv[1]
    inp = ("[input]\\nsystem=\\n   O 0.0 0.0 0.0\\n   H 0.77259794 0.55567785 0.0\\n"
           "   H -0.7731277 0.55567785 0.0\\ncharge=0\\nruntype=soc\\nbasis=6-31g\\n"
           "method=tdhf\\nfunctional=" + functional + "\\nispher=false\\nsoc_2e=1\\n\\n"
           "[scf]\\ntype=rohf\\nmultiplicity=3\\nconv=1e-10\\nmaxit=300\\n\\n"
           "[tdhf]\\ntype=mrsf\\nnstate=3\\nmultiplicity=3\\nconv=1e-9\\n")
    name = "h2o_soc"
    open(name + ".inp", "w").write(inp)
    r = Runner(project=name, input_file=name + ".inp", log=name + ".log", silent=1, usempi=False)
    r.run()
    d = r.mol.data
    nbf = int(round(math.sqrt(np.asarray(d["OQP::VEC_MO_A"]).size)))
    na = int(np.asarray(d["nelec_A"]).ravel()[0]); nb = int(np.asarray(d["nelec_B"]).ravel()[0])
    dim = na * (nbf - nb)
    out = {"nbf": nbf, "nocca": na, "noccb": nb,
           "xs": np.asarray(d["OQP::td_bvec_mo_s"]).ravel().reshape(-1, dim).tolist(),
           "xt": np.asarray(d["OQP::td_bvec_mo_t"]).ravel().reshape(-1, dim).tolist(),
           "l1e": np.asarray(d["OQP::soc_lmo_1e"]).ravel().reshape(3, nbf, nbf).transpose(0, 2, 1).tolist(),
           "l2e": np.asarray(d["OQP::soc_lmo_2e"]).ravel().reshape(3, nbf, nbf).transpose(0, 2, 1).tolist(),
           "hre": np.asarray(d["OQP::soc_hsoc_re"]).ravel().tolist(),
           "him": np.asarray(d["OQP::soc_hsoc_im"]).ravel().tolist(),
           "eval": np.asarray(d["OQP::soc_eval"]).ravel().tolist()}
    print("SOCTDM_RESULT " + json.dumps(out))
''')


def _oqp_available():
    try:
        from oqp import lib  # noqa: F401
        return hasattr(lib, "soc_mrsf")
    except Exception:
        return False


@unittest.skipUnless(_oqp_available(), "OpenQP shared library is not built")
class SOCTransitionDensities(unittest.TestCase):

    def _run(self, functional):
        with tempfile.TemporaryDirectory(prefix="soctdm_") as wd:
            drv = os.path.join(wd, "driver.py")
            with open(drv, "w") as f:
                f.write(DRIVER)
            env = dict(os.environ)
            env.setdefault("OMP_NUM_THREADS", "4")
            proc = subprocess.run([sys.executable, drv, functional], cwd=wd, env=env,
                                  capture_output=True, text=True, timeout=3600)
            lines = [l for l in proc.stdout.splitlines() if l.startswith("SOCTDM_RESULT ")]
            if proc.returncode != 0 or not lines:
                self.fail(f"driver failed:\n{proc.stdout[-3000:]}\n{proc.stderr[-3000:]}")
            return json.loads(lines[-1][len("SOCTDM_RESULT "):])

    def _check(self, functional):
        r = self._run(functional)
        xs, xt = np.array(r["xs"]), np.array(r["xt"])
        ns, nt = xs.shape[0], xt.shape[0]
        n = ns + 3 * nt
        lmo = np.array(r["l1e"]) + np.array(r["l2e"])
        hn = (np.array(r["hre"]).reshape(n, n) + 1j * np.array(r["him"]).reshape(n, n)).T
        ho = oracle_hsoc(xs, xt, r["nbf"], r["nocca"], r["noccb"], lmo)
        diff = np.abs(hn - ho).max() * DFAC
        self.assertLess(diff, 1e-8, f"native vs oracle H_SOC differ by {diff:.3e} cm-1")
        self.assertGreater(np.abs(ho).max() * DFAC, 1.0)          # a non-trivial coupling exists
        self.assertLess(np.abs(hn - hn.conj().T).max() * DFAC, 1e-10)
        # the oracle must see the former defects: a lower-core CO state couples
        self.assertGreater(np.abs(ho[:ns, ns:]).max() * DFAC, 10.0)

    def test_mrsf_hf(self):
        self._check("")

    def test_mrsf_bhhlyp(self):
        self._check("bhhlyp")


if __name__ == "__main__":
    unittest.main()
