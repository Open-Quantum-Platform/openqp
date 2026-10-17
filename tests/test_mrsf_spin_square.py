"""<S^2> printed for MRSF-TDDFT response states.

MRSF spin-adapts every configuration class except core->virtual (C->V).  The
singlet C->V pair (E+ - E-)/sqrt(2) has four open shells (C, O1, O2, V) and is
sqrt(2/3) S=0 + sqrt(1/3) S=2, so a singlet carries <S^2> = 2 w_CV while every
triplet is exactly 2 (Dier, Dittmer, Rikus, Dreuw, Neugebauer, ChemRxiv 2026,
doi:10.26434/chemrxiv.15010307, Eqs. 24-26).  OpenQP used to print a
hard-coded 0 for singlets, and the UMRSF expression kept an E+/E- coupling for
C->V pairs that does not exist, giving 1 (singlet) and 3 (triplet) per unit
C->V weight in the ROHF limit.

The reference here is a determinant-level oracle: the MRSF M_S=0 state is
built in second quantization from the packed response vector (the amplitude
map of JCP 149, 104101 and its SI), R- = S-^2 R+ / 2 is the M_S=-1 partner of
the reference triplet, and <S^2> = |S+ Psi|^2 is evaluated exactly.  The
oracle is first checked on random vectors: everything except C->V must come out
spin pure, which fixes its phase conventions independently of OpenQP.
"""
import re
import tempfile
from pathlib import Path

import pytest

try:
    import oqp                      # noqa: F401  (must precede numpy: ILP64)
    from oqp.pyoqp import Runner
    HAVE_OQP = True
except Exception:                   # pragma: no cover - source-only checkout
    HAVE_OQP = False

import numpy as np

ISQ2 = 1.0 / np.sqrt(2.0)


# ---------------------------------------------------------------------------
# Determinant oracle.  Spin orbital 2p (alpha) / 2p+1 (beta) of spatial MO p;
# a determinant is an occupation bit mask, a state a {mask: coefficient} dict.
# ---------------------------------------------------------------------------

def _op(create, mode, state):
    out = {}
    for bits, c in state.items():
        if bool((bits >> mode) & 1) == create:
            continue
        sign = -1.0 if bin(bits & ((1 << mode) - 1)).count("1") % 2 else 1.0
        key = bits ^ (1 << mode)
        out[key] = out.get(key, 0.0) + sign * c
    return out


def _add(a, b, f=1.0):
    out = dict(a)
    for k, v in b.items():
        out[k] = out.get(k, 0.0) + f * v
    return out


def _ex(state, p, q):
    """a+_p a_q"""
    return _op(True, p, _op(False, q, state))


def _s_plus(state, nmo):
    out = {}
    for p in range(nmo):
        out = _add(out, _ex(state, 2 * p, 2 * p + 1))
    return out


def _s_minus(state, nmo):
    out = {}
    for p in range(nmo):
        out = _add(out, _ex(state, 2 * p + 1, 2 * p))
    return out


def mrsf_state(x, mult, nocca, noccb, nbf):
    """MRSF M_S=0 state of the packed vector x (alpha-occupied index fastest).

    Psi = 2**-1/2 sum_{(i,a) not OO} X(i,a) [E+(i,a) + lam E-(i,a)]
        + 2**-1/2 X(O1,O1) [L + lam R]  (+ singlet: X(O2,O1) G + X(O1,O2) D)
    with E+ = a+(a beta) a(i alpha) R+, E- = a+(a alpha) a(i beta) R-,
    L = E+(O1,O1), R = E+(O2,O2), G = E+(O2,O1), D = E+(O1,O2), lam = -1/+1.
    """
    nc, o1, o2 = noccb, noccb, noccb + 1
    lam = -1.0 if mult == 1 else 1.0
    x2 = np.asarray(x, dtype=float).reshape(nbf - noccb, nocca)
    ref = {0: 1.0}
    for k in range(nc):
        ref = _op(True, 2 * k + 1, _op(True, 2 * k, ref))
    rp = _op(True, 2 * o2, _op(True, 2 * o1, ref))
    rm = {k: 0.5 * v for k, v in _s_minus(_s_minus(rp, nbf), nbf).items()}
    psi = {}
    for i in range(nocca):
        for a in range(nc, nbf):
            amp = x2[a - nc, i]
            if amp == 0.0:
                continue
            if i >= nc and a < nocca:
                if (i, a) == (o1, o1):
                    pair = _add(_ex(rp, 2 * o1 + 1, 2 * o1), _ex(rp, 2 * o2 + 1, 2 * o2), lam)
                    psi = _add(psi, pair, amp * ISQ2)
                elif mult == 1 and (i, a) == (o2, o1):
                    psi = _add(psi, _ex(rp, 2 * o1 + 1, 2 * o2), amp)
                elif mult == 1 and (i, a) == (o1, o2):
                    psi = _add(psi, _ex(rp, 2 * o2 + 1, 2 * o1), amp)
                continue
            pair = _add(_ex(rp, 2 * a + 1, 2 * i), _ex(rm, 2 * a, 2 * i + 1), lam)
            psi = _add(psi, pair, amp * ISQ2)
    return psi


def spin_square(psi, nmo):
    """<S^2> of an M_S=0 state: S^2 = S- S+ there."""
    norm = sum(v * v for v in psi.values())
    return sum(v * v for v in _s_plus(psi, nmo).values()) / norm


def _apply_s2(psi, nmo):
    """S^2 |psi> = (S- S+ + Sz^2 + Sz) |psi> for an Sz eigenstate."""
    if not psi:
        return {}
    ms = _ms(next(iter(psi)))
    out = _s_minus(_s_plus(psi, nmo), nmo)
    return _add(out, psi, ms * ms + ms)


def _ms(bits):
    n = bin(bits).count("1")
    nb = bin(bits & int("10" * (bits.bit_length() // 2 + 1), 2)).count("1")
    return 0.5 * ((n - nb) - nb)


def _dot(a, b):
    return sum(v * b.get(k, 0.0) for k, v in a.items())


def spin_weights(psi, nmo):
    """(S=0, S=1, S=2) weights of an M_S=0 state with at most four open shells:
    exact Loewdin projectors on S(S+1) in {0, 2, 6}."""
    norm = _dot(psi, psi)
    s2 = _apply_s2(psi, nmo)
    m1 = _dot(psi, s2) / norm          # <S^2>
    m2 = _dot(s2, s2) / norm           # <S^4>
    w2 = (m2 - 2.0 * m1) / 24.0        # S^2 (S^2 - 2) / 24
    w1 = (6.0 * m1 - m2) / 8.0         # S^2 (6 - S^2) / 8
    return 1.0 - w1 - w2, w1, w2


def quintet_state(x, nocca, noccb, nbf):
    """M_S=2 state of the mrst=5 vector: beta(core i) -> alpha(virtual a) on R+,
    X(i,a) packed with the core index fastest (noccb x (nbf-nocca))."""
    x2 = np.asarray(x, dtype=float).reshape(nbf - nocca, noccb)
    ref = {0: 1.0}
    for k in range(noccb):
        ref = _op(True, 2 * k + 1, _op(True, 2 * k, ref))
    rp = _op(True, 2 * (noccb + 1), _op(True, 2 * noccb, ref))
    psi = {}
    for a in range(nocca, nbf):
        for i in range(noccb):
            if x2[a - nocca, i] != 0.0:
                psi = _add(psi, _ex(rp, 2 * a, 2 * i + 1), x2[a - nocca, i])
    return psi


def spin_square_general(psi, nmo):
    """<S^2> = <S- S+> + Sz^2 + Sz for an Sz eigenstate of any M_S."""
    return _dot(psi, _apply_s2(psi, nmo)) / _dot(psi, psi)


def cv_weight(x, nocca, noccb, nbf):
    x2 = np.asarray(x, dtype=float).reshape(nbf - noccb, nocca)
    return float(np.sum(x2[2:, :noccb] ** 2) / np.sum(x2 ** 2))


def _random_vector(rng, mult, nocca, noccb, nbf, classes=("OO", "CO", "OV", "CV")):
    x2 = rng.standard_normal((nbf - noccb, nocca))
    keep = np.zeros_like(x2, dtype=bool)
    if "OO" in classes:
        keep[0:2, noccb:] = True
    if "CO" in classes:
        keep[0:2, :noccb] = True
    if "OV" in classes:
        keep[2:, noccb:] = True
    if "CV" in classes:
        keep[2:, :noccb] = True
    x2[~keep] = 0.0
    x2[1, noccb + 1] = 0.0                      # O2->O2 is not an MRSF amplitude
    if mult == 3:
        x2[0, noccb + 1] = x2[1, noccb] = 0.0   # G and D exist only for singlets
    return (x2 / np.linalg.norm(x2)).ravel()


@pytest.mark.parametrize("nocca,nbf", [(3, 6), (4, 7), (5, 8)])
def test_oracle_spin_adapted_classes_are_pure(nocca, nbf):
    rng = np.random.default_rng(nocca * 31 + nbf)
    noccb = nocca - 2
    for mult, expected in ((1, 0.0), (3, 2.0)):
        x = _random_vector(rng, mult, nocca, noccb, nbf, classes=("OO", "CO", "OV"))
        assert spin_square(mrsf_state(x, mult, nocca, noccb, nbf), nbf) == pytest.approx(
            expected, abs=1e-12)


@pytest.mark.parametrize("nocca,nbf", [(3, 6), (4, 7), (5, 8)])
def test_oracle_core_virtual_contamination(nocca, nbf):
    rng = np.random.default_rng(nocca * 17 + nbf)
    noccb = nocca - 2
    for mult in (1, 3):
        x = _random_vector(rng, mult, nocca, noccb, nbf)
        s2 = spin_square(mrsf_state(x, mult, nocca, noccb, nbf), nbf)
        expected = 2.0 * cv_weight(x, nocca, noccb, nbf) if mult == 1 else 2.0
        assert s2 == pytest.approx(expected, abs=1e-12)
        if mult == 1:
            assert s2 > 1e-3


@pytest.mark.parametrize("nocca,nbf", [(3, 6), (4, 7), (5, 8)])
def test_oracle_spin_decomposition(nocca, nbf):
    """Singlet C->V pair = sqrt(2/3) S=0 + sqrt(1/3) S=2: quintet weight w_CV/3,
    no triplet admixture; the triplet target is pure S=1."""
    rng = np.random.default_rng(nocca * 7 + nbf)
    noccb = nocca - 2
    x = _random_vector(rng, 1, nocca, noccb, nbf)
    w0, w1, w2 = spin_weights(mrsf_state(x, 1, nocca, noccb, nbf), nbf)
    wcv = cv_weight(x, nocca, noccb, nbf)
    assert (w0, w1, w2) == pytest.approx((1.0 - wcv / 3.0, 0.0, wcv / 3.0), abs=1e-12)
    x = _random_vector(rng, 3, nocca, noccb, nbf)
    assert spin_weights(mrsf_state(x, 3, nocca, noccb, nbf), nbf) == pytest.approx(
        (0.0, 1.0, 0.0), abs=1e-12)


@pytest.mark.parametrize("nocca,nbf", [(3, 6), (4, 7), (5, 8)])
def test_oracle_quintet_target_is_pure(nocca, nbf):
    """mrst=5 flips a core beta into an alpha virtual of R+: four parallel alpha
    spins, M_S = 2 = S, so <S^2> = 6 exactly for any amplitudes."""
    rng = np.random.default_rng(nocca * 13 + nbf)
    noccb = nocca - 2
    x = rng.standard_normal(noccb * (nbf - nocca))
    psi = quintet_state(x, nocca, noccb, nbf)
    assert _ms(next(iter(psi))) == 2.0
    assert spin_square_general(psi, nbf) == pytest.approx(6.0, abs=1e-12)


@pytest.mark.skipif(not HAVE_OQP, reason="compiled oqp package not importable")
def test_native_umrsf_rohf_limit_and_core_virtual_coupling():
    error = oqp.ffi.new("double *")
    checks = oqp.ffi.new("int *")
    # A missing entry point must fail instead of skipping a stale build.
    oqp.lib.mrsf_spin_square_selftest(error, checks)
    assert checks[0] == 15
    assert np.isfinite(error[0]) and error[0] < 1e-12


# ---------------------------------------------------------------------------
# End to end: CH2 MRSF/BHHLYP, printed <S^2> against the oracle
# ---------------------------------------------------------------------------

CH2_INPUT = """[input]
system=
   6   0.000000000   0.000000000   0.104000000
   1   0.000000000   0.862000000  -0.468000000
   1   0.000000000  -0.862000000  -0.468000000
charge=0
runtype=energy
basis=6-31g
functional=bhhlyp
method=tdhf

[guess]
type=huckel

[scf]
type=rohf
multiplicity=3
conv=1.0e-10
maxit=200

[tdhf]
type=mrsf
multiplicity={mult}
nstate=4
conv=1.0e-10
maxit=100
"""

S2_LINE = re.compile(r"^\s*<S\^2> =\s+(-?\d+\.\d+)\s*$", re.M)
NOCCA, NOCCB = 5, 3        # CH2 triplet ROHF reference: 8 electrons


def _state_vectors(bvec, nstate):
    """Response vectors of OQP::td_bvec_mo, one row per state.

    The Python view reports the Fortran (xvec_dim, nstates) array with shape
    (xvec_dim, nstates) but C strides over the Fortran memory, so the states
    are the rows of the C-order reshape.  Unit norm guards the layout."""
    vectors = np.asarray(bvec, dtype=float).reshape(nstate, -1)
    np.testing.assert_allclose(np.einsum("ij,ij->i", vectors, vectors), 1.0, atol=1e-8)
    return vectors


@pytest.mark.skipif(not HAVE_OQP, reason="compiled oqp package not importable")
@pytest.mark.parametrize("mult", [1, 3])
def test_ch2_mrsf_printed_spin_square_matches_determinant_oracle(mult):
    with tempfile.TemporaryDirectory(prefix="mrsf_s2_") as tmp:
        path = Path(tmp) / "ch2.inp"
        path.write_text(CH2_INPUT.format(mult=mult))
        log = path.with_suffix(".log")
        runner = Runner(input_file=str(path), log=str(log))
        runner.run()
        bvec = np.array(runner.mol.data["OQP::td_bvec_mo"], copy=True)
        printed = [float(v) for v in S2_LINE.findall(log.read_text())]

    nstate = 4
    vectors = _state_vectors(bvec, nstate)
    nbf = vectors.shape[1] // NOCCA + NOCCB
    assert NOCCA * (nbf - NOCCB) == vectors.shape[1]
    assert len(printed) == nstate, printed

    oracle = [spin_square(mrsf_state(v, mult, NOCCA, NOCCB, nbf), nbf) for v in vectors]
    np.testing.assert_allclose(printed, oracle, atol=1.5e-4)
    if mult == 1:
        # Nonzero C->V weight: the former hard-coded 0 cannot pass.
        assert max(oracle) > 1e-2, oracle
    else:
        np.testing.assert_allclose(oracle, 2.0, atol=1e-10)


@pytest.mark.skipif(not HAVE_OQP, reason="compiled oqp package not importable")
def test_ch2_quintet_printed_spin_square_matches_determinant_oracle():
    nstate = 3
    deck = CH2_INPUT.format(mult=5).replace("nstate=4", f"nstate={nstate}")
    with tempfile.TemporaryDirectory(prefix="mrsf_s2_q_") as tmp:
        path = Path(tmp) / "ch2.inp"
        path.write_text(deck)
        log = path.with_suffix(".log")
        runner = Runner(input_file=str(path), log=str(log))
        runner.run()
        bvec = np.array(runner.mol.data["OQP::td_bvec_mo"], copy=True)
        printed = [float(v) for v in S2_LINE.findall(log.read_text())]

    vectors = _state_vectors(bvec, nstate)
    nbf = vectors.shape[1] // NOCCB + NOCCA
    assert NOCCB * (nbf - NOCCA) == vectors.shape[1]
    assert len(printed) == nstate, printed
    oracle = [spin_square_general(quintet_state(v, NOCCA, NOCCB, nbf), nbf) for v in vectors]
    np.testing.assert_allclose(oracle, 6.0, atol=1e-10)
    np.testing.assert_allclose(printed, oracle, atol=1.5e-4)
