module tdhf_hessian_fixed_density_mod

  use precision, only: dp

  implicit none

  private
  public :: build_tdhf_fixed_density_hessian

contains

!###############################################################################

  subroutine build_tdhf_fixed_density_hessian(infos, h_fixed, h_ground_fixed)
    ! Fixed-density part of the total closed-shell TDHF Cartesian Hessian.
    !
    ! This is the direct second derivative of the stationary state-specific
    ! gradient with all relaxed densities held fixed.  It contains nuclear
    ! repulsion, one-electron and Pulay terms, the effective-core-potential
    ! term when present, and the two-electron integral contraction.  It does
    ! not contain the first-order orbital/amplitude/Z-vector response rows.
    ! Exchange-correlation quadrature derivatives are also excluded here.

    use types, only: information
    use basis_tools, only: basis_set
    use oqp_tagarray_driver, only: tagarray_get_data, data_has_tags, &
      OQP_DM_A, OQP_VEC_MO_A, OQP_WAO, OQP_TD_P, OQP_TD_XPY, OQP_TD_XMY
    use mathlib, only: unpack_matrix, symmetrize_matrix, orthogonal_transform
    use tdhf_lib, only: iatogen
    use tdhf_gradient_mod, only: grd2_tdhf_compute_data_t
    use grd1, only: eijden, hess_nn, hess_ee_overlap, hess_ee_kinetic, hess_en
    use grd2, only: grd2_hess_driver
    use ecp_tool, only: add_ecphess
    use messages, only: show_message, WITH_ABORT

    implicit none

    type(information), target, intent(inout) :: infos
    real(kind=dp), intent(out) :: h_fixed(:,:), h_ground_fixed(:,:)

    character(len=*), parameter :: subroutine_name = 'build_tdhf_fixed_density_hessian'
    character(len=*), parameter :: tags_required(*) = [character(len=80) :: &
      OQP_DM_A, OQP_VEC_MO_A, OQP_WAO, OQP_TD_P, OQP_TD_XPY, OQP_TD_XMY]

    type(basis_set), pointer :: basis
    type(grd2_tdhf_compute_data_t) :: gcomp
    real(kind=dp), contiguous, pointer :: dmat_a(:), mo_a(:,:), wao(:)
    real(kind=dp), contiguous, pointer :: td_p(:,:), xpy(:,:), xmy(:,:)
    real(kind=dp), allocatable, target :: d(:,:,:), p(:,:,:), xpy_ao(:,:,:), xmy_ao(:,:,:)
    real(kind=dp), allocatable :: wrk1(:,:), wrk2(:,:), overlap_density(:), one_density(:)
    integer :: natom, ncart, nbf, nocc
    real(kind=dp) :: scale_exch

    if (infos%control%scftype /= 1 .or. infos%mol_prop%mult /= 1 .or. &
        infos%tddft%mult /= 1 .or. infos%tddft%tda) then
      call show_message('TDHF fixed-density Hessian currently requires a '// &
        'closed-shell restricted reference, a singlet target, and full response.', WITH_ABORT)
    end if
    ! The XC quadrature part is assembled separately by
    ! build_tdhf_xc_fixed_hessian; this routine supplies the one- and
    ! two-electron skeleton with hybrid exchange scaling.

    basis => infos%basis
    basis%atoms => infos%atoms
    natom = size(basis%atoms%xyz, 2)
    ncart = 3*natom
    nbf = basis%nbf
    nocc = infos%mol_prop%nocc
    scale_exch = 1.0_dp
    if (infos%control%hamilton >= 20) scale_exch = infos%dft%hfscale

    if (any(shape(h_fixed) /= [ncart, ncart]) .or. &
        any(shape(h_ground_fixed) /= [ncart, ncart])) then
      call show_message('TDHF fixed-density Hessian output has the wrong Cartesian shape.', &
                        WITH_ABORT)
    end if

    call data_has_tags(infos%dat, tags_required, 'tdhf_hessian_fixed_density_mod', &
                       subroutine_name, WITH_ABORT)
    call tagarray_get_data(infos%dat, OQP_DM_A, dmat_a)
    call tagarray_get_data(infos%dat, OQP_VEC_MO_A, mo_a)
    call tagarray_get_data(infos%dat, OQP_WAO, wao)
    call tagarray_get_data(infos%dat, OQP_TD_P, td_p)
    call tagarray_get_data(infos%dat, OQP_TD_XPY, xpy)
    call tagarray_get_data(infos%dat, OQP_TD_XMY, xmy)

    allocate(d(nbf,nbf,1), p(nbf,nbf,1), xpy_ao(nbf,nbf,1), xmy_ao(nbf,nbf,1))
    allocate(wrk1(nbf,nbf), wrk2(nbf,nbf))
    allocate(overlap_density(nbf*(nbf+1)/2), one_density(nbf*(nbf+1)/2))

    call unpack_matrix(dmat_a, d(:,:,1))
    call unpack_matrix(td_p(:,1), p(:,:,1))

    call iatogen(xpy(:,infos%tddft%target_state), wrk1, nocc, nocc)
    call symmetrize_matrix(wrk1, nbf)
    wrk1 = 0.5_dp*wrk1
    call orthogonal_transform('t', nbf, mo_a, wrk1, xpy_ao(:,:,1), wrk2)

    call iatogen(xmy(:,infos%tddft%target_state), wrk1, nocc, nocc)
    call orthogonal_transform('t', nbf, mo_a, wrk1, xmy_ao(:,:,1), wrk2)

    h_fixed = 0.0_dp
    call hess_nn(basis%atoms, basis%ecp_zn_num, h_fixed)

    ! Relative to the ground-state gradient, the excited-state Lagrangian adds
    ! -2 W to the overlap density and +2 P to the one-electron density.
    call eijden(overlap_density, nbf, infos)
    overlap_density = overlap_density - 2.0_dp*wao
    one_density = dmat_a + 2.0_dp*td_p(:,1)

    call hess_ee_overlap(basis, overlap_density, h_fixed)
    call hess_ee_kinetic(basis, one_density, h_fixed)
    call hess_en(basis, basis%atoms%xyz, &
                 basis%atoms%zn - basis%ecp_zn_num, one_density, h_fixed)
    call add_ecphess(basis, basis%atoms%xyz, one_density, h_fixed)

    gcomp = grd2_tdhf_compute_data_t(d2=d, p2=p, xpy2=xpy_ao, xmy2=xmy_ao, &
                                      hfscale=scale_exch, nbf=nbf)
    call gcomp%init()
    call gcomp%build_cart(basis)
    call grd2_hess_driver(infos, basis, h_fixed, gcomp)
    call gcomp%clean()

    ! Ground-state skeleton.  The complete ground-state response is taken from
    ! the independently validated native RHF Hessian; keeping this skeleton
    ! separate lets the caller replace only its response part and avoids the
    ! noncanonical occupied-block ambiguity in dW.
    h_ground_fixed = 0.0_dp
    call hess_nn(basis%atoms, basis%ecp_zn_num, h_ground_fixed)
    call eijden(overlap_density, nbf, infos)
    one_density = dmat_a
    call hess_ee_overlap(basis, overlap_density, h_ground_fixed)
    call hess_ee_kinetic(basis, one_density, h_ground_fixed)
    call hess_en(basis, basis%atoms%xyz, &
                 basis%atoms%zn - basis%ecp_zn_num, one_density, h_ground_fixed)
    call add_ecphess(basis, basis%atoms%xyz, one_density, h_ground_fixed)
    p = 0.0_dp; xpy_ao = 0.0_dp; xmy_ao = 0.0_dp
    gcomp = grd2_tdhf_compute_data_t(d2=d, p2=p, xpy2=xpy_ao, xmy2=xmy_ao, &
                                      hfscale=scale_exch, nbf=nbf)
    call gcomp%init()
    call gcomp%build_cart(basis)
    call grd2_hess_driver(infos, basis, h_ground_fixed, gcomp)
    call gcomp%clean()

    deallocate(d, p, xpy_ao, xmy_ao, wrk1, wrk2, overlap_density, one_density)

  end subroutine build_tdhf_fixed_density_hessian

end module tdhf_hessian_fixed_density_mod
