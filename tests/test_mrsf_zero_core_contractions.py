"""Regression guard for two-electron MRSF references.

For H2, ``nocca == 2`` and the closed-shell part of ``umrsfcbc`` has zero
width.  Accelerate BLAS rejects the resulting zero-leading-dimension array
section even though the contraction is mathematically empty.
"""

from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MRSF_LIB = ROOT / "source" / "tdhf_mrsf_lib.F90"
NAC_DRIVER = ROOT / "source" / "modules" / "mrsf_nac_driver.F90"
NAC_METRIC = ROOT / "source" / "modules" / "mrsf_nac_metric_data.F90"
STATE_OVERLAP = ROOT / "source" / "modules" / "get_states_overlap.F90"
NAC_KERNEL = ROOT / "pyoqp" / "oqp" / "library" / "nac_kernel.py"
TDXC_GRAD = ROOT / "source" / "dftlib" / "dft_gridint_tdxc_grad.F90"
FXC_GRID = ROOT / "source" / "dftlib" / "dft_gridint_fxc.F90"
CPHF = ROOT / "source" / "modules" / "cphf.F90"
NAC_INTERCHANGE = ROOT / "source" / "modules" / "mrsf_nac_interchange.F90"


def test_umrsfcbc_skips_empty_closed_shell_blas_contraction():
    source = MRSF_LIB.read_text()
    start = source.index("subroutine umrsfcbc(")
    end = source.index("end subroutine umrsfcbc", start)
    body = source[start:end]

    contraction = body.index("call dgemm('n','t', nbf, nocca-2")
    guard = body.rindex("if (nocca > 2) then", 0, contraction)
    close = body.index("end if", contraction)

    assert guard < contraction < close
    assert "zero-width array" in body[max(0, guard - 240):contraction]


def test_mrsfmntoia_skips_zero_leading_dimension_outputs():
    source = MRSF_LIB.read_text()
    start = source.index("subroutine mrsfmntoia(")
    end = source.index("end subroutine mrsfmntoia", start)
    body = source[start:end]

    calls = "call dgemm('t','n',noca-2,1,nbf"
    assert body.count(calls) == 2
    assert body.count("if (noca > 2) then") == 2
    for section in body.split(calls)[1:]:
        assert "end if" in section[:240]


def test_resident_nac_accepts_a_zero_closed_shell_core():
    source = NAC_DRIVER.read_text()
    assert "nocb64 < 0_c_int64_t" in source
    assert "0 <= nocb < noca <= nbf and two SOMOs" in source
    assert "nocb64 < 1_c_int64_t" not in source

    metric = NAC_METRIC.read_text()
    assert "noca - nocb /= 2 .or. nocb < 0" in metric
    assert "with a core" not in metric
    assert "at least one doubly occupied orbital" not in metric


def test_coreless_sia_minor_is_the_literal_one_by_one_limit():
    """The final valid ov_exact case(3) write defines M(2,a) for H2."""
    metric = NAC_METRIC.read_text()
    start = metric.index("pure subroutine build_ia_maps(")
    end = metric.index("end subroutine build_ia_maps", start)
    maps = metric[start:end]
    assert "if (noc == 1) then" in maps
    assert "rows(1) = noc + 1" in maps
    assert "cols(1) = nocb + j1" in maps
    assert maps.index("if (noc == 1) then") < maps.index("do k = 1, noc - 2")

    overlap = STATE_OVERLAP.read_text()
    case3 = overlap.index("case (3)", overlap.index("subroutine ov_exact("))
    end = overlap.index("end select", case3)
    body = overlap[case3:end]
    assert "if (noc == 1) then" in body
    assert "temp1 = s_mo(noc+1,ia1)" in body
    assert body.index("if (noc == 1) then") < body.index("!  (1,1) block")

    kernel = NAC_KERNEL.read_text()
    start = kernel.index("def _ia_maps(")
    end = kernel.index("def s_ia_of(", start)
    maps_py = kernel[start:end]
    assert "if noc == 1:" in maps_py
    assert "return [noc], [nocb + j1]" in maps_py


def test_disabled_fxc_does_not_alias_unallocated_gradient_scratch():
    source = TDXC_GRAD.read_text()
    start = source.index("subroutine update(self, xce, mythread)")
    end = source.index("end subroutine", start)
    update = source[start:end]
    associate = update[update.index("associate ("):update.index(")", update.index("associate ("))]
    assert "grad_x =>" not in associate
    assert update.count("self%grad_x(:,:,1,j,mythread)") == 1
    assert update.count("self%grad_x(:,:,2,j,mythread)") == 1

    parallel = source[source.index("subroutine parallel_start("):start]
    assert "self%rtau(nSpin, xce%maxPts, self%nMtx, nthreads)" in parallel
    assert "if (xce%funTyp == OQP_FUNTYP_MGGA) then\n        allocate" not in parallel

    fxc = FXC_GRID.read_text()
    parallel = fxc[fxc.index("subroutine parallel_start("):fxc.index("end subroutine", fxc.index("subroutine parallel_start("))]
    assert "self%drRho(3, nSpin, xce%maxPts, self%nMtx, nThreads)" in parallel
    assert "self%drRho(4," not in parallel
    assert "self%rTau(nSpin, xce%maxPts, self%nMtx, nThreads)" in parallel
    mgga = parallel[parallel.index("if (xce%funTyp == OQP_FUNTYP_MGGA) then"):]
    assert "self%rTau(" not in mgga
    assert "self%moG1_(" in mgga


def test_rohf_cphf_skips_zero_extent_beta_blocks():
    source = CPHF.read_text()
    scalar = source[source.index("subroutine cphf_apbx_rohf("):
                    source.index("end subroutine cphf_apbx_rohf", source.index("subroutine cphf_apbx_rohf("))]
    batch = source[source.index("subroutine cphf_apbx_rohf_batch("):
                   source.index("end subroutine cphf_apbx_rohf_batch", source.index("subroutine cphf_apbx_rohf_batch("))]
    assert scalar.count("if (noccb > 0) then") == 3
    assert "dm = 0.0_dp" in scalar
    assert batch.count("if (noccb > 0) then") == 3
    assert "p%dxb_batch(:,:,ivec) = 0.0_dp" in batch


def test_rohf_nac_adjoint_skips_coreless_beta_blas_blocks():
    source = NAC_INTERCHANGE.read_text()
    routines = (
        "mrsf_nac_rohf_hf_adjoint",
        "mrsf_nac_rohf_hf_adjoint_batch",
        "mrsf_nac_xc_adjoint",
        "mrsf_nac_xc_adjoint_batch",
    )
    for name in routines:
        start = source.index(f"subroutine {name}(")
        end = source.index(f"end subroutine {name}", start)
        body = source[start:end]
        beta_density = body.index("pzb", body.index("if (noccb > 0) then"))
        assert "if (noccb > 0) then" in body[:beta_density]
        assert "pzb" in body[beta_density:]
        assert "= 0.0_dp" in body[beta_density:body.index("end if", beta_density) + 6]

    batch_start = source.index("subroutine mrsf_nac_xc_adjoint_batch(")
    helper = source.index("subroutine ao_to_mo_occ(", batch_start)
    helper_end = source.index("end subroutine ao_to_mo_occ", helper)
    helper_body = source[helper:helper_end]
    assert "if (nocc == 0) return" in helper_body
    assert helper_body.index("if (nocc == 0) return") < helper_body.index("call dgemm")
