"""Nonorthogonal UMRSF S² against an independent determinant expansion.

The production code contracts overlap cofactors. This reference applies creation,
annihilation and a physical spin rotation to determinants in a common spatial
orthonormal basis. It contains no overlap cofactor/SVD formula.
"""
import os
from pathlib import Path
import subprocess
import sys

import pytest
from test_mrsf_spin_square import (
    _add, _op, _random_vector, spin_square, CH2_INPUT, S2_LINE, PRINT_ATOL,
    _state_vectors, oqp, OQP_IMPORT_ERROR, spin_square_general,
)
import numpy as np


def _orbital_op(state, vector, spin, create):
    out = {}
    for p, coefficient in enumerate(vector):
        if coefficient != 0:
            out = _add(out, _op(create, 2*p+spin, state), coefficient)
    return out


def _spin_flip(state):
    out = {}
    for bits, coefficient in state.items():
        term = {0: coefficient}
        for p in reversed(range(bits.bit_length())):
            if bits >> p & 1:
                term = _op(True, p ^ 1, term)
                if p % 2:
                    term = {d: -c for d, c in term.items()}
        out = _add(out, term)
    return out


def nonorthogonal_state(a, b, nc, x, mult):
    """Explicit paired-parent extension; unrestricted parents have spin parity,
    not necessarily pure S. R-=Ry(pi)R+, never S-²R+/2 for a contaminated R+."""
    n = len(a); na = nc+2
    rp = {0: 1.0}
    for i in range(nc):
        rp = _orbital_op(_orbital_op(rp, a[:, i], 0, True), b[:, i], 1, True)
    for i in (nc, nc+1):
        rp = _orbital_op(rp, a[:, i], 0, True)
    rm = _spin_flip(rp)
    lam = -1 if mult == 1 else 1

    def pair(i, j):
        plus = _orbital_op(_orbital_op(rp, a[:, i], 0, False), b[:, j], 1, True)
        minus = _orbital_op(_orbital_op(rm, a[:, i], 1, False), b[:, j], 0, True)
        return _add(plus, minus, lam)

    psi = {}
    amplitudes = np.asarray(x).reshape(n-nc, na)
    for j in range(nc, n):
        for i in range(na):
            value = amplitudes[j-nc, i]
            if value == 0:
                continue
            if i >= nc and j < na:
                if (i, j) == (nc, nc):
                    vec = _add(pair(nc, nc), pair(nc+1, nc+1), lam)
                    value /= 2*np.sqrt(2)
                elif mult == 1 and i != j:
                    vec = pair(i, j)
                    value /= 2
                else:
                    continue
            else:
                vec = pair(i, j)
                value /= np.sqrt(2)
            psi = _add(psi, vec, value)
    return psi


def nonorthogonal_quintet(a, b, nc, x):
    """Explicit Ms=2 beta-core -> alpha-virtual determinant expansion."""
    n = len(a); na = nc+2
    rp = {0: 1.0}
    for i in range(nc):
        rp = _orbital_op(_orbital_op(rp, a[:, i], 0, True), b[:, i], 1, True)
    for i in (nc, nc+1):
        rp = _orbital_op(rp, a[:, i], 0, True)
    amplitudes = np.asarray(x).reshape(n-na, nc)
    psi = {}
    for j in range(na, n):
        for i in range(nc):
            term = _orbital_op(_orbital_op(rp, b[:, i], 1, False), a[:, j], 0, True)
            psi = _add(psi, term, amplitudes[j-na, i])
    return psi


def native_s2(a, b, nc, x, mult, metric=None):
    assert oqp is not None, f"compiled oqp unavailable: {OQP_IMPORT_ERROR!r}"
    n = len(a)
    if metric is None:
        metric = np.eye(n)
    packed_s = np.array([metric[i, j] for j in range(n) for i in range(j+1)])
    arrays = [np.array(y, dtype=np.float64, order='F').ravel(order='F')
              for y in (a, b, packed_s, x)]
    out = oqp.ffi.new('double *')
    oqp.lib.mrsf_spin_square_evaluate(n, nc, mult,
        *(oqp.ffi.from_buffer('double[]', y) for y in arrays), out)
    return out[0]


@pytest.mark.parametrize('theta', [0.0, 0.2, 0.6, 1.0, np.pi/2])
@pytest.mark.parametrize('mult', [1, 3])
def test_added_beta_virtual_overlap(theta, mult):
    # B_C2 wedge B_V2 = C2 wedge V2. The old diagonal formula gives
    # 2+sin(theta)^2 instead of 2 for C1->V2, in BOTH multiplicities.
    a = np.eye(6); b = a.copy()
    b[:, [1, 5]] = b[:, [1, 5]] @ [[np.cos(theta), -np.sin(theta)],
                                  [np.sin(theta), np.cos(theta)]]
    x = np.zeros(16); x[12] = 1
    assert spin_square(nonorthogonal_state(a, b, 2, x, mult), 6) == pytest.approx(2, abs=1e-12)
    assert native_s2(a, b, 2, x, mult) == pytest.approx(2, abs=1e-12)


@pytest.mark.parametrize('nc,n', [(0, 2), (1, 3), (2, 4), (0, 4), (1, 5), (2, 6), (3, 7)])
@pytest.mark.parametrize('mult', [1, 3])
def test_nonorthogonal_classes_and_cross_terms(nc, n, mult):
    rng = np.random.default_rng(515+nc)
    a = np.linalg.qr(rng.normal(size=(n, n)))[0]
    b = np.linalg.qr(rng.normal(size=(n, n)))[0]
    na = nc+2
    mixed = _random_vector(rng, mult, na, nc, n)
    # Every occupied slot and every coupling to a nonsymmetric mixed vector.
    for index in np.flatnonzero(mixed):
        x = np.zeros_like(mixed); x[index] = 1
        for trial in (x, mixed+x):
            ref = spin_square(nonorthogonal_state(a, b, nc, trial, mult), n)
            assert native_s2(a, b, nc, trial, mult) == pytest.approx(ref, abs=1e-12)
    for scale in (1e-150, -3.7, 1e150):
        ref = spin_square(nonorthogonal_state(a, b, nc, mixed, mult), n)
        assert native_s2(a, b, nc, mixed*scale, mult) == pytest.approx(ref, abs=1e-12)


@pytest.mark.parametrize('rank', [0, 1, 2])
@pytest.mark.parametrize('mult', [1, 3])
def test_zero_singular_values_and_degenerate_complement(rank, mult):
    n, nc = 6, 2
    rng = np.random.default_rng(710+rank)
    a = np.eye(n); b = np.eye(n)
    # Exact zeros (permutations), not a cos(pi/2) roundoff surrogate.
    for j in range(nc-rank):
        b[:, [j, 4+j]] = b[:, [4+j, j]]
    x = _random_vector(rng, mult, nc+2, nc, n)
    ref = spin_square(nonorthogonal_state(a, b, nc, x, mult), n)
    assert native_s2(a, b, nc, x, mult) == pytest.approx(ref, abs=1e-12)
    # Common arbitrary spatial coordinate rotation must leave the result alone.
    q = np.linalg.qr(rng.normal(size=(n, n)))[0]
    assert native_s2(q@a, q@b, nc, x, mult) == pytest.approx(ref, abs=1e-12)


@pytest.mark.parametrize('mult', [1, 3])
def test_nonidentity_ao_metric(mult):
    rng = np.random.default_rng(1515)
    n, nc = 6, 2
    a = np.linalg.qr(rng.normal(size=(n, n)))[0]
    b = np.linalg.qr(rng.normal(size=(n, n)))[0]
    transform = rng.normal(size=(n, n))*0.1+np.diag(np.arange(1, n+1))
    metric = transform.T@transform
    x = _random_vector(rng, mult, nc+2, nc, n)
    ref = spin_square(nonorthogonal_state(a, b, nc, x, mult), n)
    assert native_s2(np.linalg.solve(transform, a), np.linalg.solve(transform, b),
                     nc, x, mult, metric) == pytest.approx(ref, abs=1e-12)


@pytest.mark.parametrize('case', ['zero', 'nonfinite', 'null', 'near_null', 'spin'])
def test_invalid_paired_state_fails(case):
    assert oqp is not None, f"compiled oqp unavailable: {OQP_IMPORT_ERROR!r}"
    code = '''
from test_umrsf_spin_square_nonorthogonal import native_s2
import numpy as np
n=3; a=np.eye(n); b=a.copy(); x=np.zeros(6); mult=3
x[0]=1
if CASE == 'zero': x[:]=0
elif CASE == 'nonfinite': x[0]=np.nan
elif CASE in ('null','near_null'):
    b[:,[0,1]]=b[:,[1,0]]
    if CASE == 'near_null': x[4]=1e-9
elif CASE == 'spin': mult=7
print(native_s2(a,b,0,x,mult))
'''.replace('CASE', repr(case))
    env = dict(os.environ)
    env['PYTHONPATH'] = str(Path(__file__).parent)
    run = subprocess.run([sys.executable, '-c', code], env=env, text=True, capture_output=True)
    assert run.returncode != 0
    assert 'UMRSF spin square:' in run.stdout+run.stderr


@pytest.mark.parametrize('molecule,nc', [('ch2', 3), ('lih', 1)])
@pytest.mark.parametrize('mult', [1, 3, 5])
def test_umrsf_printed_states_against_common_basis_oracle(tmp_path, molecule, nc, mult):
    from oqp.pyoqp import Runner
    deck = (CH2_INPUT.format(mult=mult).replace('basis=6-31g', 'basis=sto-3g')
            .replace('type=rohf', 'type=uhf').replace('type=mrsf', 'type=umrsf'))
    if molecule == 'lih':
        start = deck.index('system='); end = deck.index('charge=')
        deck = deck[:start]+'system=\n  Li 0 0 0\n  H 0 0 3.0\n'+deck[end:]
    nstate = 3 if mult == 5 else 4
    deck = deck.replace('nstate=4', f'nstate={nstate}')
    path = tmp_path / f'{molecule}.inp'; path.write_text(deck)
    log = path.with_suffix('.log')
    runner = Runner(input_file=str(path), log=str(log))
    runner.run()
    data = runner.mol.data
    a = np.array(data['OQP::VEC_MO_A'], copy=True).T
    b = np.array(data['OQP::VEC_MO_B'], copy=True).T
    n = len(a)
    metric = np.zeros((n, n)); metric[np.tril_indices(n)] = np.asarray(data['OQP::SM']).ravel()
    metric += np.tril(metric, -1).T
    # Cholesky is just an AO coordinate map here: L^T C has Euclidean metric.
    # No change to the orbital or response spaces is made.
    transform = np.linalg.cholesky(metric).T
    aa, bb = transform@a, transform@b
    np.testing.assert_allclose(aa.T@aa, np.eye(n), atol=1e-10)
    np.testing.assert_allclose(bb.T@bb, np.eye(n), atol=1e-10)
    vectors = _state_vectors(np.array(data['OQP::td_bvec_mo'], copy=True), nstate)
    printed = np.array([float(v) for v in S2_LINE.findall(log.read_text())])
    assert len(printed) == nstate
    assert np.max(np.abs(aa.T@bb-np.eye(n))) > 1e-3
    if mult == 5:
        exact = [spin_square_general(nonorthogonal_quintet(aa, bb, nc, x), n) for x in vectors]
    else:
        exact = [spin_square(nonorthogonal_state(aa, bb, nc, x, mult), n) for x in vectors]
    native = [native_s2(a, b, nc, x, mult, metric) for x in vectors]
    np.testing.assert_allclose(native, exact, rtol=0, atol=1e-10)
    np.testing.assert_allclose(printed, exact, rtol=0, atol=PRINT_ATOL)


@pytest.mark.parametrize('rank', [0, 1, 2])
@pytest.mark.parametrize('mult', [1, 3])
def test_occupied_and_core_rotations_change_svd_gauge(rank, mult):
    rng = np.random.default_rng(902+rank)
    n, nc, na = 6, 2, 4
    a = np.eye(n); b = np.eye(n)
    for j in range(nc-rank):
        b[:, [j, 4+j]] = b[:, [4+j, j]]
    # Only virtual particles: arbitrary occupied-alpha rotations need no OO fold.
    amplitudes = np.zeros((na, n-nc)); amplitudes[:, 2:] = rng.normal(size=(na, 2))
    x = amplitudes.ravel(order='F')
    ref = spin_square(nonorthogonal_state(a, b, nc, x, mult), n)
    occupied_rotation = np.linalg.qr(rng.normal(size=(na, na)))[0]
    core_rotation = np.linalg.qr(rng.normal(size=(nc, nc)))[0]
    a[:, :na] = a[:, :na]@occupied_rotation
    b[:, :nc] = b[:, :nc]@core_rotation
    x = (occupied_rotation.T@amplitudes).ravel(order='F')
    rotated_ref = spin_square(nonorthogonal_state(a, b, nc, x, mult), n)
    assert rotated_ref == pytest.approx(ref, abs=1e-12)
    assert native_s2(a, b, nc, x, mult) == pytest.approx(ref, abs=1e-12)


@pytest.mark.parametrize('theta', [0.0, 0.2, 0.6, 1.0, np.pi/2])
def test_quintet_contaminated_parent(theta):
    a = np.eye(6); b = a.copy()
    b[:, [1, 5]] = b[:, [1, 5]] @ [[np.cos(theta), -np.sin(theta)],
                                  [np.sin(theta), np.cos(theta)]]
    x = np.array([1.0, 0.0, 0.0, 0.0])  # beta C1 -> alpha V1
    ref = 6+np.sin(theta)**2
    assert spin_square_general(nonorthogonal_quintet(a, b, 2, x), 6) == pytest.approx(ref, abs=1e-12)
    assert native_s2(a, b, 2, x, 5) == pytest.approx(ref, abs=1e-12)


@pytest.mark.parametrize('nc,n', [(1, 4), (1, 5), (2, 5), (2, 6), (3, 7)])
def test_quintet_cross_terms_metric_and_rotations(nc, n):
    rng = np.random.default_rng(2515+nc+n)
    na = nc+2; nv = n-na
    a = np.linalg.qr(rng.normal(size=(n, n)))[0]
    b = np.linalg.qr(rng.normal(size=(n, n)))[0]
    mixed = rng.normal(size=nc*nv)
    for index in range(nc*nv):
        unit = np.zeros(nc*nv); unit[index] = 1
        for x in (unit, mixed+unit):
            ref = spin_square_general(nonorthogonal_quintet(a, b, nc, x), n)
            assert native_s2(a, b, nc, x, 5) == pytest.approx(ref, abs=1e-12)
            if nc == 1:
                assert ref == pytest.approx(6, abs=1e-12)
    ref = spin_square_general(nonorthogonal_quintet(a, b, nc, mixed), n)
    for scale in (1e-150, -3.7, 1e150):
        assert native_s2(a, b, nc, mixed*scale, 5) == pytest.approx(ref, abs=1e-12)
    transform = rng.normal(size=(n, n))*0.1+np.diag(np.arange(1, n+1))
    assert native_s2(np.linalg.solve(transform, a), np.linalg.solve(transform, b),
                     nc, mixed, 5, transform.T@transform) == pytest.approx(ref, abs=1e-12)
    # Rotate both excitation indices and counter-transform amplitudes. The
    # occupied-alpha determinant changes only by an irrelevant overall phase.
    u = np.linalg.qr(rng.normal(size=(nc, nc)))[0]
    v = np.linalg.qr(rng.normal(size=(nv, nv)))[0]
    occupied = np.linalg.qr(rng.normal(size=(na, na)))[0]
    a[:, :na] = a[:, :na]@occupied
    a[:, na:] = a[:, na:]@v
    b[:, :nc] = b[:, :nc]@u
    rotated = (u.T@mixed.reshape(nv, nc).T@v).ravel(order='F')
    assert spin_square_general(nonorthogonal_quintet(a, b, nc, rotated), n) == pytest.approx(ref, abs=1e-12)
    assert native_s2(a, b, nc, rotated, 5) == pytest.approx(ref, abs=1e-12)
    # Restricted limit, including a nonsymmetric mixed vector.
    assert native_s2(a, a, nc, mixed, 5) == pytest.approx(6, abs=1e-12)


@pytest.mark.parametrize('case', ['zero', 'nonfinite', 'no_core', 'no_virtual'])
def test_invalid_quintet_state_fails(case):
    assert oqp is not None, f"compiled oqp unavailable: {OQP_IMPORT_ERROR!r}"
    code = """
from test_umrsf_spin_square_nonorthogonal import native_s2
import numpy as np
n=6; nc=2; x=np.ones(4)
if CASE == 'zero': x[:]=0
elif CASE == 'nonfinite': x[0]=np.nan
elif CASE == 'no_core': nc=0
elif CASE == 'no_virtual': n=4
print(native_s2(np.eye(n),np.eye(n),nc,x,5))
""".replace('CASE', repr(case))
    env = dict(os.environ)
    env['PYTHONPATH'] = str(Path(__file__).parent)
    run = subprocess.run([sys.executable, '-c', code], env=env, text=True, capture_output=True)
    assert run.returncode != 0
    assert 'UMRSF spin square:' in run.stdout+run.stderr
