"""Static contracts for the single-call resident MRSF NAC driver."""

from pathlib import Path
import sys
from types import SimpleNamespace

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
DRIVER = ROOT / "source" / "modules" / "mrsf_nac_driver.F90"
METRIC = ROOT / "source" / "modules" / "mrsf_nac_metric_data.F90"
INTERCHANGE = ROOT / "source" / "modules" / "mrsf_nac_interchange.F90"
INT2 = ROOT / "source" / "integrals" / "int2.F90"
SCF_ADDONS = ROOT / "source" / "scf_addons.F90"
FOCK_DERIV = ROOT / "source" / "modules" / "fock_deriv.F90"
GRD2 = ROOT / "source" / "integrals" / "grd2.F90"
GRD2_RYS = ROOT / "source" / "integrals" / "grd2_rys.F90"
HEADER = ROOT / "include" / "oqp.h"
PYTHON = ROOT / "pyoqp" / "oqp" / "library" / "nac_analytic.py"


def test_single_production_entry_is_exported_and_python_is_thin():
    driver = DRIVER.read_text()
    header = HEADER.read_text()
    production = PYTHON.read_text()
    assert 'bind(C, name="mrsf_nac_lagrangian")' in driver
    assert "void mrsf_nac_lagrangian(struct oqp_handle_t *inf);" in header
    assert production.count("oqp.mrsf_nac_lagrangian(mol)") == 1
    for forbidden in (
        "oqp.mrsf_nac_metric_data(mol)",
        "oqp.mrsf_nac_wpair(mol",
        "oqp.mrsf_nac_amp_pair(mol",
        "oqp.mrsf_nac_response(mol)",
        "oqp.mrsf_nac_rohf_zvector(mol)",
        "for I in range(nstate)",
        "for J in range(nstate)",
        "NAC_ANALYTIC_FORWARD_GATE",
        "oqp.hf_hessian(mol)",
    ):
        assert forbidden not in production


def test_driver_uses_actual_resident_state_count_and_streamed_metric():
    driver = DRIVER.read_text()
    metric = METRIC.read_text()
    assert "nstate64 = infos%tddft%nstate" in driver
    assert "size(energies) /= nstate" in driver
    assert "size(bvec_mo,2) /= nstate" in driver
    assert "call mrsf_nac_metric_column(infos, jstate, gamma_column)" in driver
    assert "call mrsf_nac_metric_data" not in driver
    assert "raw_sij(noca,noca,nstate)" in metric
    assert "raw_sab(nvirb,nvirb,nstate)" in metric
    assert "raw_sia(noca,nvirb,nstate)" in metric
    assert "energies_saved = energies" in driver
    assert "gap = energies_saved(pair_j(batch_pair))" in driver
    assert "energies_saved(pair_i(batch_pair))" in driver
    assert "default_int_limit64/state_pair_size64" in driver


def test_driver_batches_one_adjoint_per_unordered_pair():
    driver = DRIVER.read_text()
    body = driver.split(
        "subroutine mrsf_nac_lagrangian(infos, gradient_rhs, gradient_solution)", 1
    )[1].split("end subroutine mrsf_nac_lagrangian", 1)[0]
    ordered_source_sequence = (
        "call mrsf_nac_wpair_batch_impl(",
        "mt_frozen_tag = reshape(wpair_mt(:,:,wpair_index)",
        "call mrsf_nac_amp(infos, istate, jstate)",
        "call mrsf_nac_esum(infos, istate, jstate)",
        "call mrsf_nac_response(infos)",
        "call mrsf_nac_rohf_pair_overlap(infos)",
    )
    positions = [body.index(token) for token in ordered_source_sequence]
    assert positions == sorted(positions)
    assert "npair = nstate*(nstate - 1)/2" in body
    assert "pair_sign = merge(0.5_dp, -0.5_dp, istate < jstate)" in body
    assert "gamma_pair = pair_sign*gamma_column(:,istate)" in body
    assert "rhs_batch(:,ipair) = rhs_batch(:,ipair) + rhs_in" in body
    assert "pair_sign*rhs_in" not in body
    assert "nonz_batch(coord,ipair) = nonz_batch(coord,ipair) +" in body
    assert body.count("call mrsf_nac_wpair_batch_impl(") == 1
    assert "integer, parameter :: wpair_batch_width = 3" in body
    assert "call mrsf_nac_wpair_impl(infos, istate, jstate)" not in body
    assert body.count("call mrsf_nac_amp(infos, istate, jstate)") == 1
    assert body.count("call mrsf_nac_esum(infos, istate, jstate)") == 1
    assert body.count("call mrsf_nac_response(infos)") == 1
    assert body.count("call mrsf_nac_rohf_pair_overlap(infos)") == 1
    assert body.count(
        "call mrsf_nac_rohf_pair_overlap(infos, metric_only=.true.)"
    ) == 1
    direct_guard = body.index("if (istate < jstate) then")
    direct_call = body.index("call mrsf_nac_wpair_batch_impl(")
    metric_publish = body.index("gamma_tag = gamma_pair", direct_call)
    assert direct_guard < direct_call < metric_publish
    assert body.count("call mrsf_nac_rohf_zvector_batch(") == 1
    # One exact contraction is always present; the second is guarded by the
    # research-only predictor audit and cannot replace the production vector.
    assert body.count("call mrsf_nac_xc_adjoint_batch(") == 2
    assert body.count("call mrsf_nac_rohf_hf_adjoint_batch(") == 2
    assert "if (audit_enabled .and. any(predictor_available)) then" in body
    assert "integer, parameter :: z_batch_width = 3" in body
    assert "do z_first = 1, npair, z_batch_width" in body
    assert "rhs_batch(:,z_first:z_last)" in body
    assert "solution_batch(:,z_first:z_last)" in body
    assert "integer, parameter :: hf_batch_width = 3" in body
    assert "integer, parameter :: xc_batch_width = 3" in body
    assert "do xc_first = 1, npair, xc_batch_width" in body
    assert "solution_batch(:,xc_first:xc_last)" in body
    assert "xc_batch(:,:,xc_first:xc_last)" in body
    assert "solution_batch(:,hf_first:hf_last)" in body
    assert "hf_batch(:,:,hf_first:hf_last)" in body
    assert "call mrsf_nac_rohf_hf_adjoint(infos)" not in body
    assert "call mrsf_nac_xc_adjoint(infos)" not in body
    assert body.count("call mrsf_nac_pair_accumulate_antisym(") == 1

    pair_loop = body.index("do istate = 1, nstate")
    pair_skip = body.index("if (istate == jstate) cycle", pair_loop)
    pair_source = body.index("call mrsf_nac_rohf_pair_overlap(infos)")
    batch_call = body.index("call mrsf_nac_rohf_zvector_batch(")
    hf_batch_call = body.index("call mrsf_nac_rohf_hf_adjoint_batch(")
    xc_batch_call = body.index("call mrsf_nac_xc_adjoint_batch(")
    adjoint_loop = body.index("do ipair = 1, npair", batch_call)
    hf_publish = body.index("hf_tag = hf_batch(:,:,ipair)")
    xc_publish = body.index("xc_tag = xc_batch(:,:,ipair)")
    accumulate_call = body.index("call mrsf_nac_pair_accumulate_antisym(")
    assert (
        pair_loop < pair_skip < pair_source < batch_call < hf_batch_call
        < xc_batch_call < adjoint_loop < hf_publish < xc_publish
        < accumulate_call
    )
    assert "call mrsf_nac_rohf_zvector(infos)" not in body
    assert body.count("call cphf_solve_rohf") == 1
    assert body.index("if (present(gradient_rhs)) then") < body.index(
        "call cphf_solve_rohf"
    )
    assert "hf_hessian" not in body


def test_hf_adjoint_batch_shares_pair_independent_work_and_jk_pass():
    source = INTERCHANGE.read_text()
    body = source.split(
        "subroutine mrsf_nac_rohf_hf_adjoint_batch(infos, z_vectors, ghf_vectors)",
        1,
    )[1].split("end subroutine mrsf_nac_rohf_hf_adjoint_batch", 1)[0]
    assert body.count("call der_overlap_matrix(basis, dsa)") == 1
    assert body.count("call der_kinetic_matrix(basis, dta)") == 1
    assert body.count("call der_nucattr_matrix(") == 1
    assert body.count("call ecp_deriv_ints(") == 1
    assert body.count("call fock_jk(") == 1
    assert body.count("call fock_deriv_contract_os_batch(") == 1
    assert "call fock_deriv_contract_os(" not in body
    assert "allocate(dmz(nbf2,2*nrhs), vjkz(nbf2,2*nrhs)" in body
    assert "do irhs = 1, nrhs" in body

    int2 = INT2.read_text()
    scf_addons = SCF_ADDONS.read_text()
    update = int2.split(
        "subroutine int2_urohf_data_t_update(this, buf)", 1
    )[1].split("end subroutine int2_urohf_data_t_update", 1)[0]
    assert "this%nfocks = size(this%d,2)" in int2
    assert "do ifock = 1, this%nfocks, 2" in update
    assert "sum(this%d(ij,ifock:ifock+1))" in update
    assert scf_addons.count("int2_urohf_data_t(nfocks=size(d,2)") == 2


def test_hf_derivative_eri_batch_shares_recurrence_without_nested_openmp():
    fock = FOCK_DERIV.read_text()
    scalar = fock.split(
        "subroutine fock_deriv_contract_os(infos, basis, pcoul, pexch, mmat, hfscale, gx)",
        1,
    )[1].split("end subroutine fock_deriv_contract_os", 1)[0]
    assert scalar.count("call grd2_driver(infos, basis, gx, gcomp)") == 1
    assert "grd2_driver_batch" not in scalar
    fock_batch = fock.split(
        "subroutine fock_deriv_contract_os_batch(", 1
    )[1].split("end subroutine fock_deriv_contract_os_batch", 1)[0]
    assert fock_batch.count("call grd2_driver_batch(") == 1
    assert "integer, parameter :: max_rhs = 3" in fock_batch
    assert "gcomps(ia)%pexch => pexcha" in fock_batch
    assert "gcomps(ib)%pexch => pexchb" in fock_batch
    assert "gall(:,:,2*irhs-1) + gall(:,:,2*irhs)" in fock_batch

    grd2 = GRD2.read_text()
    driver = grd2.split(
        "subroutine grd2_driver_batch_gen(infos, basis, de, gcomps)", 1
    )[1].split("end subroutine grd2_driver_batch_gen", 1)[0]
    assert driver.count("!$omp parallel") == 1
    assert "reduction(+:skip1, skip2, numint)" in driver
    assert "reduction(+:skip1, skip2, numint, de)" not in driver
    assert "call batch_worker(skip1, skip2, numint)" in driver
    assert "recursive subroutine batch_worker(" in driver
    parallel_region = driver.split("!$omp parallel", 1)[1].split(
        "!$omp end parallel", 1
    )[0]
    assert "allocatable" not in parallel_region
    assert "private(" not in parallel_region
    worker = driver.split("subroutine batch_worker(", 1)[1].split(
        "end subroutine batch_worker", 1
    )[0]
    assert "allocatable :: dab(:,:), dabmax(:), fd_batch(:,:,:)" in worker
    assert "allocatable :: de_thread(:,:,:)" in worker
    assert "type(grd2_int_data_t) :: gdat" in worker
    assert "!$omp do schedule(dynamic,4) collapse(2)" in worker
    assert "de_thread(:,gdat%at,iprobe)" in worker
    assert "!$omp critical(grd2_batch_de_merge)" in worker
    assert "de = de + de_thread" in worker
    assert "probe_active = dabmax*gmax*real(q4,dp) >= cutoff2" in worker
    assert "nquartet = product(basis%naos(gdat%id))" in worker
    assert "dab(1:nquartet,iprobe) = 0.0_dp" in worker
    assert "product(gdat%nbf)" not in worker
    assert driver.count("call grd2_rys_compute_batch(") == 2
    assert "do iprobe = 2, nprobe" in driver
    assert "gcomps(iprobe)%attenuated .neqv. gcomps(1)%attenuated" in driver
    assert "gcomps(iprobe)%mu /= gcomps(1)%mu" in driver
    assert "uniform attenuation and mu" in driver

    rys = GRD2_RYS.read_text()
    recurrence = rys.split(
        "subroutine compute_grd_ints_batch(", 1
    )[1].split("end subroutine compute_grd_ints_batch", 1)[0]
    assert recurrence.count("call compute_rys_rw(") == 1
    assert recurrence.count("call compute_coefficients(") == 1
    assert recurrence.count("call compute_der_xyz_ijkl(") == 1
    assert recurrence.count("call compute_der_ijkl_batch(") == 1
    assert "!$omp parallel" not in recurrence


def test_python_production_call_cannot_enable_forward_cphf(monkeypatch):
    import importlib.util

    spec = importlib.util.spec_from_file_location(
        "nac_analytic_no_forward_cphf", PYTHON
    )
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)

    class Data(dict):
        def __init__(self):
            super().__init__({
                "OQP::td_bvec_mo": np.arange(6.0).reshape(2, 3),
                "OQP::td_energies": np.array([0.1, 0.2, 0.3]),
                "natom": 1,
            })
            self._data = SimpleNamespace(
                control=SimpleNamespace(int2e_cutoff=1.0e-12),
                tddft=SimpleNamespace(nstate=3),
            )

    mol = SimpleNamespace(
        config={
            "scf": {"conv": 1.0e-10},
            "tdhf": {"conv": 1.0e-10, "multiplicity": 1},
        },
        data=Data(),
    )
    calls = []

    def resident_driver(active_mol):
        calls.append("mrsf_nac_lagrangian")
        active_mol.data["OQP::nac_dcv"] = np.zeros(27)
        active_mol.data["OQP::nac_nacv"] = np.zeros(27)

    def forbidden_forward_driver(_):
        raise AssertionError("production analytic NAC called forward CPHF")

    fake_oqp = SimpleNamespace(
        mrsf_nac_lagrangian=resident_driver,
        hf_hessian=forbidden_forward_driver,
    )
    monkeypatch.setitem(sys.modules, "oqp", fake_oqp)
    monkeypatch.setenv("NAC_ANALYTIC_FORWARD_GATE", "1")
    nacv, dcv = module.analytic_nac(mol)

    assert calls == ["mrsf_nac_lagrangian"]
    assert nacv.shape == (3, 3, 1, 3)
    assert dcv.shape == (3, 3, 1, 3)


def test_driver_guards_scope_gaps_and_restores_mutated_state():
    driver = DRIVER.read_text()
    for required in (
        "infos%control%scftype /= 3",
        "infos%tddft%umrsf",
        "infos%tddft%mult /= 1",
        "infos%control%conv > 1.0e-8_dp",
        "infos%tddft%cnvtol > 1.0e-8_dp",
        ".not. infos%mol_energy%SCF_converged",
        ".not. infos%mol_energy%Davidson_converged",
        "abs(gap) <= gap_floor",
        "bvec_mo = bvec_saved",
        "infos%control%int2e_cutoff = cutoff_saved",
    ):
        assert required in driver
    assert "inquire(unit=iw, opened=log_was_open)" in driver
