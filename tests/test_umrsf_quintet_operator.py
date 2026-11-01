"""Beta->alpha quintet response against a complete independent AO-integral matrix."""
import oqp  # native package is required, before NumPy for the ILP64 runtime
import numpy as np
import pytest

DECK = """[input]
system=
 C 0 0 0.104
 H 0 0.862 -0.468
 H 0 -0.862 -0.468
basis={basis}
runtype=energy
method=tdhf
functional={functional}
[guess]
type=huckel
[scf]
type={reference}
multiplicity=3
conv=1e-10
maxit=200
[tdhf]
type={response}
multiplicity=5
nstate={states}
conv=1e-10
maxit=100
"""


def _unpack(data, n):
    matrix = np.zeros((n, n))
    matrix[np.tril_indices(n)] = np.asarray(data).ravel()
    return matrix+np.tril(matrix, -1).T


def quintet_matrix(a, b, fa, fb, eri, nc, exchange):
    """A_ia,jb=delta_ij Fa_ab-delta_ab Fb_ij-cHF (Bi Bj|Aa Ab)."""
    na = nc+2; nv = len(a)-na
    faa = a[:, na:].T@fa@a[:, na:]
    fbb = b[:, :nc].T@fb@b[:, :nc]
    two = np.einsum('pi,qj,ra,sb,pqrs->aibj', b[:, :nc], b[:, :nc],
                    a[:, na:], a[:, na:], eri, optimize=True)
    result = np.einsum('ij,ab->aibj', np.eye(nc), faa)
    result -= np.einsum('ab,ij->aibj', np.eye(nv), fbb)
    result -= exchange*two
    return result.reshape(nc*nv, nc*nv)


@pytest.mark.parametrize('reference,response,functional,exchange,basis,dimension', [
    ('uhf', 'umrsf', '', 1.0, 'sto-3g', 6),
    ('rohf', 'mrsf', '', 1.0, 'sto-3g', 6),
    ('uhf', 'umrsf', 'bhhlyp', 0.5, 'sto-3g', 6),
    ('rohf', 'mrsf', 'bhhlyp', 0.5, 'sto-3g', 6),
    ('uhf', 'umrsf', 'bhhlyp', 0.5, 'cc-pvdz', 57),
])
def test_complete_quintet_operator(tmp_path, reference, response, functional, exchange, basis, dimension):
    # Native availability is required; a missing installed package fails.
    from oqp.pyoqp import Runner
    from oqp.library.ints_2e import eri_ao
    path = tmp_path/'ch2.inp'
    path.write_text(DECK.format(reference=reference, response=response, functional=functional, basis=basis, states=dimension))
    runner = Runner(input_file=str(path), log=str(path.with_suffix('.log')))
    runner.run()
    assert runner.mol.mol_energy.Davidson_converged
    data = runner.mol.data
    a = np.array(data['OQP::VEC_MO_A'], copy=True).T
    b = np.array(data['OQP::VEC_MO_B'], copy=True).T
    assert data['nelec_A'] == 5 and data['nelec_B'] == 3
    n = len(a); nc = int(data['nelec_B'])
    assert dimension == nc*(n-nc-2)
    fa = _unpack(data['OQP::FOCK_A'], n)
    fb = _unpack(data['OQP::FOCK_B'], n)
    eri = eri_ao(runner.mol)
    for array in (a, b, fa, fb, eri):
        assert np.isfinite(array).all()
    expected = quintet_matrix(a, b, fa, fb, eri, nc, exchange)
    assert np.isfinite(expected).all()
    np.testing.assert_allclose(expected, expected.T, rtol=0, atol=1e-12, equal_nan=False)
    values = np.array(data['OQP::td_energies'], copy=True).ravel()
    vectors = np.array(data['OQP::td_bvec_mo'], copy=True).reshape(dimension, dimension).T
    assert np.isfinite(values).all() and np.isfinite(vectors).all()
    np.testing.assert_allclose(vectors.T@vectors, np.eye(dimension), rtol=0, atol=1e-10, equal_nan=False)
    np.testing.assert_allclose(values, np.linalg.eigvalsh(expected), rtol=0, atol=1e-9, equal_nan=False)
    np.testing.assert_allclose(expected@vectors, vectors*values, rtol=0, atol=1e-9, equal_nan=False)
    # Complete set of Ritz vectors recovers every signed matrix element.
    np.testing.assert_allclose((vectors*values)@vectors.T, expected, rtol=0, atol=1e-9, equal_nan=False)
    rng = np.random.default_rng(515)
    uc = np.linalg.qr(rng.normal(size=(nc, nc)))[0]
    uv = np.linalg.qr(rng.normal(size=(n-nc-2, n-nc-2)))[0]
    ar = a.copy(); br = b.copy()
    for orbitals in (ar, br):
        orbitals[:, :nc] = orbitals[:, :nc]@uc
        orbitals[:, nc+2:] = orbitals[:, nc+2:]@uv
    rotated = quintet_matrix(ar, br, fa, fb, eri, nc, exchange)
    transform = np.kron(uv, uc)
    assert np.isfinite(rotated).all() and np.isfinite(transform).all()
    np.testing.assert_allclose(rotated, transform.T@expected@transform, rtol=0, atol=1e-12, equal_nan=False)
    # Exercise the production driver in the rotated basis, not only the oracle.
    # Both reference densities stay unchanged, so the AO Fock and ERIs remain valid.
    # Stored orbital energies enter the initial vectors/preconditioner, but the
    # complete response space is spanned here and sigma uses the full MO Fock.
    data['OQP::VEC_MO_A'][:] = ar.T
    data['OQP::VEC_MO_B'][:] = br.T
    getattr(oqp, f'tdhf_{response}_energy')(runner.mol)
    assert runner.mol.mol_energy.Davidson_converged
    rotated_values = np.array(data['OQP::td_energies'], copy=True).ravel()
    rotated_vectors = np.array(data['OQP::td_bvec_mo'], copy=True).reshape(dimension, dimension).T
    assert np.isfinite(rotated_values).all() and np.isfinite(rotated_vectors).all()
    np.testing.assert_allclose(rotated_vectors.T@rotated_vectors, np.eye(dimension),
                               rtol=0, atol=1e-10, equal_nan=False)
    np.testing.assert_allclose(rotated_values, values, rtol=0, atol=1e-9, equal_nan=False)
    np.testing.assert_allclose(rotated@rotated_vectors, rotated_vectors*rotated_values,
                               rtol=0, atol=1e-9, equal_nan=False)
    np.testing.assert_allclose((rotated_vectors*rotated_values)@rotated_vectors.T, rotated,
                               rtol=0, atol=1e-9, equal_nan=False)
    if reference == 'uhf':
        assert np.max(abs(a-b)) > 1e-3
    else:
        np.testing.assert_allclose(a, b, rtol=0, atol=1e-12, equal_nan=False)


@pytest.mark.parametrize('reference,response', [('uhf', 'umrsf'), ('rohf', 'mrsf')])
def test_cam_quintet_rotation_covariance(tmp_path, reference, response):
    """Native CAM operator covariance, separate from the bare-ERI oracle above.

    CAM exchange is alpha/r + beta*erf(mu*r)/r; eri_ao exports only 1/r.
    Reconstruct the complete native matrix, then compare its transformation
    with an independent driver run on rotated core and virtual orbitals.
    """
    from oqp.pyoqp import Runner
    path = tmp_path/'ch2-cam.inp'
    path.write_text(DECK.format(reference=reference, response=response,
                               functional='cam-b3lyp', basis='sto-3g', states=6))
    runner = Runner(input_file=str(path), log=str(path.with_suffix('.log')))
    runner.run()
    data = runner.mol.data
    assert data['nelec_A'] == 5 and data['nelec_B'] == 3
    assert data['hamilton'] == 20
    assert data['dft'].cam_flag
    np.testing.assert_allclose([data['cam_alpha'], data['cam_beta'], data['cam_mu']],
                               [0.19, 0.46, 0.33], rtol=0, atol=1e-12, equal_nan=False)

    def eigensystem():
        assert runner.mol.mol_energy.Davidson_converged
        e = np.array(data['OQP::td_energies'], copy=True).ravel()
        x = np.array(data['OQP::td_bvec_mo'], copy=True).reshape(6, 6).T
        assert e.shape == (6,)
        assert np.isfinite(e).all() and np.isfinite(x).all()
        np.testing.assert_allclose(x.T@x, np.eye(6), rtol=0, atol=1e-10, equal_nan=False)
        return e, x

    e, x = eigensystem()
    a = np.array(data['OQP::VEC_MO_A'], copy=True).T
    b = np.array(data['OQP::VEC_MO_B'], copy=True).T
    assert np.isfinite(a).all() and np.isfinite(b).all()
    if reference == 'uhf':
        assert np.max(abs(a-b)) > 1e-3
    else:
        np.testing.assert_allclose(a, b, rtol=0, atol=1e-12, equal_nan=False)
    rng = np.random.default_rng(522)
    uc = np.linalg.qr(rng.normal(size=(3, 3)))[0]
    uv = np.linalg.qr(rng.normal(size=(2, 2)))[0]
    transform = np.kron(uv, uc)
    expected = transform.T@((x*e)@x.T)@transform
    assert np.isfinite(expected).all()
    # AO densities/Fock are invariant. Stored energies only seed/precondition
    # the solver; all six response directions are present in this test.
    for key, orbitals in [('OQP::VEC_MO_A', a), ('OQP::VEC_MO_B', b)]:
        orbitals[:, :3] = orbitals[:, :3]@uc
        orbitals[:, 5:] = orbitals[:, 5:]@uv
        data[key][:] = orbitals.T
    getattr(oqp, f'tdhf_{response}_energy')(runner.mol)
    er, xr = eigensystem()
    np.testing.assert_allclose(er, e, rtol=0, atol=1e-9, equal_nan=False)
    np.testing.assert_allclose(expected@xr, xr*er, rtol=0, atol=1e-9, equal_nan=False)
    np.testing.assert_allclose((xr*er)@xr.T, expected, rtol=0, atol=1e-9, equal_nan=False)
