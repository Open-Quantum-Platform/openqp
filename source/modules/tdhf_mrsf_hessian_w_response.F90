module tdhf_mrsf_hessian_w_response_mod

  use precision, only: dp

  implicit none
  private

  integer, parameter, public :: MRSF_W_SUCCESS = 0
  integer, parameter, public :: MRSF_W_BAD_DIMENSION = -1
  integer, parameter, public :: MRSF_W_BAD_REFERENCE = -2
  integer, parameter, public :: MRSF_W_BAD_TARGET = -3
  integer, parameter, public :: MRSF_W_UMRSF_UNSUPPORTED = -4
  integer, parameter, public :: MRSF_W_CAM_UNSUPPORTED = -5
  integer, parameter, public :: MRSF_W_ALLOCATION_FAILURE = -6

  public :: build_mrsf_w_ao
  public :: build_mrsf_w_ao_derivative
  public :: mrsf_w_response_status_message

contains

  ! Hiroya Nakata's TDHF/TDDFT analytical Hessian formulation is the
  ! methodological starting point for this MRSF extension.
  !
  ! This routine reproduces the spin-adapted two-SOMO mrsfrowcal mapping.
  ! Its result includes both factors of one half applied by the ROHF caller:
  !
  !   W_AO = (1/4) [ C W_MO C^T + (C W_MO C^T)^T ].
  !
  ! Singlet and triplet target dependence is already contained in the supplied
  ! spin-adapted response intermediates.  No additional state transformation is
  ! introduced here.
  subroutine build_mrsf_w_ao(mo, mo_energy, fa, fb, z, hxa, hxb, ppija, ppijb, &
                             noca, nocb, rohf_reference, reference_multiplicity, &
                             target_multiplicity, umrsf, cam_flag, w_ao, status)

    real(kind=dp), intent(in) :: mo(:,:), mo_energy(:), fa(:,:), fb(:,:)
    real(kind=dp), intent(in) :: z(:), hxa(:,:), hxb(:,:)
    real(kind=dp), intent(in) :: ppija(:,:), ppijb(:,:)
    integer, intent(in) :: noca, nocb, reference_multiplicity
    integer, intent(in) :: target_multiplicity
    logical, intent(in) :: rohf_reference, umrsf, cam_flag
    real(kind=dp), intent(out) :: w_ao(:,:)
    integer, intent(out) :: status

    real(kind=dp), allocatable :: scr(:,:), wrk(:,:), wmo(:,:), transformed(:,:)
    integer :: allocation_status, nbf

    w_ao = 0.0_dp
    nbf = size(mo, 1)
    if (size(w_ao, 1) /= nbf .or. size(w_ao, 2) /= nbf) then
      status = MRSF_W_BAD_DIMENSION
      return
    end if
    call validate_baseline_inputs(mo, mo_energy, fa, fb, z, hxa, hxb, ppija, &
                                  ppijb, noca, nocb, status)
    if (status /= MRSF_W_SUCCESS) return

    call validate_method(rohf_reference, reference_multiplicity, &
                         target_multiplicity, umrsf, cam_flag, status)
    if (status /= MRSF_W_SUCCESS) return

    allocate(scr(nbf,nbf), wrk(nbf,nbf), wmo(nbf,nbf), &
             transformed(nbf,nbf), stat=allocation_status)
    if (allocation_status /= 0) then
      status = MRSF_W_ALLOCATION_FAILURE
      return
    end if

    call unpack_response(z, noca, nocb, scr)
    call build_w_mo_kernel(mo_energy, fa, fb, scr, hxa, hxb, ppija, ppijb, &
                           noca, nocb, wrk, wmo)
    transformed = matmul(matmul(mo, wmo), transpose(mo))
    w_ao = 0.25_dp*(transformed + transpose(transformed))

    deallocate(scr, wrk, wmo, transformed)
  end subroutine build_mrsf_w_ao


  ! Exact first differential of the expression above.  Every product is
  ! differentiated as d(A B)=(dA)B+A(dB); products among first differentials
  ! are never formed.  The last index enumerates independent perturbations.
  subroutine build_mrsf_w_ao_derivative( &
      mo, dmo, mo_energy, dmo_energy, fa, dfa, fb, dfb, z, dz, &
      hxa, dhxa, hxb, dhxb, ppija, dppija, ppijb, dppijb, noca, nocb, &
      rohf_reference, reference_multiplicity, target_multiplicity, umrsf, &
      cam_flag, dw_ao, status)

    real(kind=dp), intent(in) :: mo(:,:), dmo(:,:,:)
    real(kind=dp), intent(in) :: mo_energy(:), dmo_energy(:,:)
    real(kind=dp), intent(in) :: fa(:,:), dfa(:,:,:), fb(:,:), dfb(:,:,:)
    real(kind=dp), intent(in) :: z(:), dz(:,:)
    real(kind=dp), intent(in) :: hxa(:,:), dhxa(:,:,:)
    real(kind=dp), intent(in) :: hxb(:,:), dhxb(:,:,:)
    real(kind=dp), intent(in) :: ppija(:,:), dppija(:,:,:)
    real(kind=dp), intent(in) :: ppijb(:,:), dppijb(:,:,:)
    integer, intent(in) :: noca, nocb, reference_multiplicity
    integer, intent(in) :: target_multiplicity
    logical, intent(in) :: rohf_reference, umrsf, cam_flag
    real(kind=dp), intent(out) :: dw_ao(:,:,:)
    integer, intent(out) :: status

    real(kind=dp), allocatable :: scr(:,:), dscr(:,:), wrk(:,:)
    real(kind=dp), allocatable :: wmo(:,:), dwmo(:,:), transformed(:,:)
    integer :: allocation_status, nbf, npert, perturbation

    dw_ao = 0.0_dp
    nbf = size(mo, 1)
    npert = size(dmo, 3)

    call validate_baseline_inputs(mo, mo_energy, fa, fb, z, hxa, hxb, ppija, &
                                  ppijb, noca, nocb, status)
    if (status /= MRSF_W_SUCCESS) return

    call validate_method(rohf_reference, reference_multiplicity, &
                         target_multiplicity, umrsf, cam_flag, status)
    if (status /= MRSF_W_SUCCESS) return

    call validate_derivative_inputs( &
        dmo, dmo_energy, dfa, dfb, dz, dhxa, dhxb, dppija, dppijb, &
        nbf, noca, nocb, npert, dw_ao, status)
    if (status /= MRSF_W_SUCCESS) return

    allocate(scr(nbf,nbf), dscr(nbf,nbf), wrk(nbf,nbf), wmo(nbf,nbf), &
             dwmo(nbf,nbf), transformed(nbf,nbf), stat=allocation_status)
    if (allocation_status /= 0) then
      status = MRSF_W_ALLOCATION_FAILURE
      return
    end if

    call unpack_response(z, noca, nocb, scr)
    call build_w_mo_kernel(mo_energy, fa, fb, scr, hxa, hxb, ppija, ppijb, &
                           noca, nocb, wrk, wmo)

    do perturbation = 1, npert
      call unpack_response(dz(:,perturbation), noca, nocb, dscr)
      call differentiate_w_mo_kernel( &
          mo_energy, dmo_energy(:,perturbation), fa, dfa(:,:,perturbation), &
          fb, dfb(:,:,perturbation), scr, dscr, dhxa(:,:,perturbation), &
          dhxb(:,:,perturbation), dppija(:,:,perturbation), &
          dppijb(:,:,perturbation), noca, nocb, wrk, dwmo)

      transformed = matmul(matmul(dmo(:,:,perturbation), wmo), transpose(mo)) &
                  + matmul(matmul(mo, dwmo), transpose(mo)) &
                  + matmul(matmul(mo, wmo), transpose(dmo(:,:,perturbation)))
      dw_ao(:,:,perturbation) = 0.25_dp*(transformed + transpose(transformed))
    end do

    deallocate(scr, dscr, wrk, wmo, dwmo, transformed)
  end subroutine build_mrsf_w_ao_derivative


  subroutine validate_baseline_inputs(mo, mo_energy, fa, fb, z, hxa, hxb, &
                                      ppija, ppijb, noca, nocb, status)

    real(kind=dp), intent(in) :: mo(:,:), mo_energy(:), fa(:,:), fb(:,:)
    real(kind=dp), intent(in) :: z(:), hxa(:,:), hxb(:,:)
    real(kind=dp), intent(in) :: ppija(:,:), ppijb(:,:)
    integer, intent(in) :: noca, nocb
    integer, intent(out) :: status

    integer :: nbf, nz

    status = MRSF_W_BAD_DIMENSION
    nbf = size(mo, 1)
    if (nbf < 1 .or. size(mo, 2) /= nbf) return
    if (nocb < 1 .or. noca - nocb /= 2 .or. noca >= nbf) return
    nz = nocb*(2 + nbf - noca) + 2*(nbf - noca)

    if (size(mo_energy) /= nbf) return
    if (size(fa, 1) /= nbf .or. size(fa, 2) /= nbf) return
    if (size(fb, 1) /= nbf .or. size(fb, 2) /= nbf) return
    if (size(z) /= nz) return
    if (size(hxa, 1) /= nbf .or. size(hxa, 2) /= noca) return
    if (size(hxb, 1) /= nbf .or. size(hxb, 2) /= nbf) return
    if (size(ppija, 1) /= noca .or. size(ppija, 2) /= noca) return
    if (size(ppijb, 1) /= nocb .or. size(ppijb, 2) /= nocb) return

    status = MRSF_W_SUCCESS
  end subroutine validate_baseline_inputs


  subroutine validate_method(rohf_reference, reference_multiplicity, &
                             target_multiplicity, umrsf, cam_flag, status)

    logical, intent(in) :: rohf_reference, umrsf, cam_flag
    integer, intent(in) :: reference_multiplicity, target_multiplicity
    integer, intent(out) :: status

    if (.not. rohf_reference .or. reference_multiplicity /= 3) then
      status = MRSF_W_BAD_REFERENCE
    else if (target_multiplicity /= 1 .and. target_multiplicity /= 3) then
      status = MRSF_W_BAD_TARGET
    else if (umrsf) then
      status = MRSF_W_UMRSF_UNSUPPORTED
    else if (cam_flag) then
      status = MRSF_W_CAM_UNSUPPORTED
    else
      status = MRSF_W_SUCCESS
    end if
  end subroutine validate_method


  subroutine validate_derivative_inputs( &
      dmo, dmo_energy, dfa, dfb, dz, dhxa, dhxb, dppija, dppijb, &
      nbf, noca, nocb, npert, output, status)

    real(kind=dp), intent(in) :: dmo(:,:,:), dmo_energy(:,:)
    real(kind=dp), intent(in) :: dfa(:,:,:), dfb(:,:,:), dz(:,:)
    real(kind=dp), intent(in) :: dhxa(:,:,:), dhxb(:,:,:)
    real(kind=dp), intent(in) :: dppija(:,:,:), dppijb(:,:,:)
    real(kind=dp), intent(in) :: output(:,:,:)
    integer, intent(in) :: nbf, noca, nocb, npert

    integer, intent(out) :: status
    integer :: nz

    status = MRSF_W_BAD_DIMENSION
    nz = nocb*(2 + nbf - noca) + 2*(nbf - noca)
    if (npert < 1) return
    if (.not. shape_is_3d(dmo, nbf, nbf, npert)) return
    if (size(dmo_energy, 1) /= nbf .or. size(dmo_energy, 2) /= npert) return
    if (.not. shape_is_3d(dfa, nbf, nbf, npert)) return
    if (.not. shape_is_3d(dfb, nbf, nbf, npert)) return
    if (size(dz, 1) /= nz .or. size(dz, 2) /= npert) return
    if (.not. shape_is_3d(dhxa, nbf, noca, npert)) return
    if (.not. shape_is_3d(dhxb, nbf, nbf, npert)) return
    if (.not. shape_is_3d(dppija, noca, noca, npert)) return
    if (.not. shape_is_3d(dppijb, nocb, nocb, npert)) return
    if (.not. shape_is_3d(output, nbf, nbf, npert)) return

    status = MRSF_W_SUCCESS
  end subroutine validate_derivative_inputs


  logical function shape_is_3d(array, first, second, third)
    real(kind=dp), intent(in) :: array(:,:,:)
    integer, intent(in) :: first, second, third

    shape_is_3d = size(array, 1) == first .and. size(array, 2) == second &
               .and. size(array, 3) == third
  end function shape_is_3d


  subroutine unpack_response(z, noca, nocb, scr)
    real(kind=dp), intent(in) :: z(:)
    integer, intent(in) :: noca, nocb
    real(kind=dp), intent(out) :: scr(:,:)

    integer :: index, i, j, a, nbf

    scr = 0.0_dp
    nbf = size(scr, 1)
    index = 0

    do i = nocb + 1, noca
      do j = 1, nocb
        index = index + 1
        scr(j,i) = z(index)
      end do
    end do

    do a = noca + 1, nbf
      do j = 1, nocb
        index = index + 1
        scr(j,a) = z(index)
      end do
    end do

    do a = noca + 1, nbf
      do i = nocb + 1, noca
        index = index + 1
        scr(i,a) = z(index)
      end do
    end do
  end subroutine unpack_response


  subroutine build_w_mo_kernel(mo_energy, fa, fb, scr, hxa, hxb, ppija, ppijb, &
                               noca, nocb, wrk, wmo)

    real(kind=dp), intent(in) :: mo_energy(:), fa(:,:), fb(:,:), scr(:,:)
    real(kind=dp), intent(in) :: hxa(:,:), hxb(:,:), ppija(:,:), ppijb(:,:)
    integer, intent(in) :: noca, nocb
    real(kind=dp), intent(inout) :: wrk(:,:)
    real(kind=dp), intent(out) :: wmo(:,:)

    integer :: a, b, i, j, k, nbf, open_first, open_last, x, y

    nbf = size(fa, 1)
    open_first = nocb + 1
    open_last = noca
    wrk = 0.0_dp
    wmo = 0.0_dp

    do x = 1, nocb
      do k = 1, nocb
        wrk(x,1:2) = wrk(x,1:2) - fa(k,x)*scr(k,open_first:open_last)
      end do
      do a = noca + 1, nbf
        wrk(x,1:2) = wrk(x,1:2) + scr(open_first:open_last,a)*fa(a,x)
      end do
    end do
    wmo(1:nocb,open_first:open_last) = 0.5_dp*wrk(1:nocb,1:2) &
        + hxa(1:nocb,open_first:open_last) &
        + hxb(1:nocb,open_first:open_last) &
        + ppija(1:nocb,open_first:open_last)
    do x = open_first, open_last
      wmo(1:nocb,x) = wmo(1:nocb,x) + mo_energy(1:nocb)*scr(1:nocb,x)
    end do

    wrk = 0.0_dp
    do i = 1, nocb
      do a = noca + 1, nbf
        wrk(i,a-noca) = fa(open_first,i)*scr(open_first,a) &
                       + fa(open_last,i)*scr(open_last,a)
      end do
    end do
    wmo(1:nocb,noca+1:nbf) = 0.5_dp*wrk(1:nocb,1:nbf-noca) &
                            + hxb(1:nocb,noca+1:nbf)
    do a = noca + 1, nbf
      wmo(1:nocb,a) = wmo(1:nocb,a) + mo_energy(1:nocb)*scr(1:nocb,a)
    end do

    wrk = 0.0_dp
    do a = noca + 1, nbf
      do k = 1, nocb
        wrk(1,a-noca) = wrk(1,a-noca) + fa(k,open_first)*scr(k,a)
        wrk(2,a-noca) = wrk(2,a-noca) + fa(k,open_last)*scr(k,a)
      end do
      wrk(1,a-noca) = wrk(1,a-noca) &
          - fb(open_first,open_first)*scr(open_first,a) &
          - fb(open_last,open_first)*scr(open_last,a)
      wrk(2,a-noca) = wrk(2,a-noca) &
          - fb(open_first,open_last)*scr(open_first,a) &
          - fb(open_last,open_last)*scr(open_last,a)
    end do
    wmo(open_first:open_last,noca+1:nbf) = 0.5_dp*wrk(1:2,1:nbf-noca) &
                                         + hxb(open_first:open_last,noca+1:nbf)
    do a = noca + 1, nbf
      wmo(open_first:open_last,a) = wmo(open_first:open_last,a) &
          + mo_energy(open_first:open_last)*scr(open_first:open_last,a)
    end do

    do i = 1, nocb
      do j = 1, i
        wmo(i,j) = ppija(i,j) + ppijb(i,j) + hxa(j,i)
      end do
    end do
    do x = open_first, open_last
      do y = open_first, x
        wmo(x,y) = hxa(y,x) + hxb(y,x) + ppija(x,y)
      end do
    end do
    do a = noca + 1, nbf
      do b = noca + 1, a
        wmo(a,b) = hxb(b,a)
      end do
    end do

    do i = 1, nbf
      wmo(i,i) = 0.5_dp*wmo(i,i)
    end do
    wmo = -wmo
  end subroutine build_w_mo_kernel


  subroutine differentiate_w_mo_kernel( &
      mo_energy, dmo_energy, fa, dfa, fb, dfb, scr, dscr, dhxa, dhxb, &
      dppija, dppijb, noca, nocb, wrk, dwmo)

    real(kind=dp), intent(in) :: mo_energy(:), dmo_energy(:)
    real(kind=dp), intent(in) :: fa(:,:), dfa(:,:), fb(:,:), dfb(:,:)
    real(kind=dp), intent(in) :: scr(:,:), dscr(:,:)
    real(kind=dp), intent(in) :: dhxa(:,:), dhxb(:,:)
    real(kind=dp), intent(in) :: dppija(:,:), dppijb(:,:)
    integer, intent(in) :: noca, nocb
    real(kind=dp), intent(inout) :: wrk(:,:)
    real(kind=dp), intent(out) :: dwmo(:,:)

    integer :: a, b, i, j, k, nbf, open_first, open_last, x, y

    nbf = size(fa, 1)
    open_first = nocb + 1
    open_last = noca
    wrk = 0.0_dp
    dwmo = 0.0_dp

    do x = 1, nocb
      do k = 1, nocb
        wrk(x,1:2) = wrk(x,1:2) &
            - dfa(k,x)*scr(k,open_first:open_last) &
            - fa(k,x)*dscr(k,open_first:open_last)
      end do
      do a = noca + 1, nbf
        wrk(x,1:2) = wrk(x,1:2) &
            + dscr(open_first:open_last,a)*fa(a,x) &
            + scr(open_first:open_last,a)*dfa(a,x)
      end do
    end do
    dwmo(1:nocb,open_first:open_last) = 0.5_dp*wrk(1:nocb,1:2) &
        + dhxa(1:nocb,open_first:open_last) &
        + dhxb(1:nocb,open_first:open_last) &
        + dppija(1:nocb,open_first:open_last)
    do x = open_first, open_last
      dwmo(1:nocb,x) = dwmo(1:nocb,x) &
          + dmo_energy(1:nocb)*scr(1:nocb,x) &
          + mo_energy(1:nocb)*dscr(1:nocb,x)
    end do

    wrk = 0.0_dp
    do i = 1, nocb
      do a = noca + 1, nbf
        wrk(i,a-noca) = dfa(open_first,i)*scr(open_first,a) &
                       + fa(open_first,i)*dscr(open_first,a) &
                       + dfa(open_last,i)*scr(open_last,a) &
                       + fa(open_last,i)*dscr(open_last,a)
      end do
    end do
    dwmo(1:nocb,noca+1:nbf) = 0.5_dp*wrk(1:nocb,1:nbf-noca) &
                             + dhxb(1:nocb,noca+1:nbf)
    do a = noca + 1, nbf
      dwmo(1:nocb,a) = dwmo(1:nocb,a) &
          + dmo_energy(1:nocb)*scr(1:nocb,a) &
          + mo_energy(1:nocb)*dscr(1:nocb,a)
    end do

    wrk = 0.0_dp
    do a = noca + 1, nbf
      do k = 1, nocb
        wrk(1,a-noca) = wrk(1,a-noca) &
            + dfa(k,open_first)*scr(k,a) + fa(k,open_first)*dscr(k,a)
        wrk(2,a-noca) = wrk(2,a-noca) &
            + dfa(k,open_last)*scr(k,a) + fa(k,open_last)*dscr(k,a)
      end do
      wrk(1,a-noca) = wrk(1,a-noca) &
          - dfb(open_first,open_first)*scr(open_first,a) &
          - fb(open_first,open_first)*dscr(open_first,a) &
          - dfb(open_last,open_first)*scr(open_last,a) &
          - fb(open_last,open_first)*dscr(open_last,a)
      wrk(2,a-noca) = wrk(2,a-noca) &
          - dfb(open_first,open_last)*scr(open_first,a) &
          - fb(open_first,open_last)*dscr(open_first,a) &
          - dfb(open_last,open_last)*scr(open_last,a) &
          - fb(open_last,open_last)*dscr(open_last,a)
    end do
    dwmo(open_first:open_last,noca+1:nbf) = 0.5_dp*wrk(1:2,1:nbf-noca) &
                                          + dhxb(open_first:open_last,noca+1:nbf)
    do a = noca + 1, nbf
      dwmo(open_first:open_last,a) = dwmo(open_first:open_last,a) &
          + dmo_energy(open_first:open_last)*scr(open_first:open_last,a) &
          + mo_energy(open_first:open_last)*dscr(open_first:open_last,a)
    end do

    do i = 1, nocb
      do j = 1, i
        dwmo(i,j) = dppija(i,j) + dppijb(i,j) + dhxa(j,i)
      end do
    end do
    do x = open_first, open_last
      do y = open_first, x
        dwmo(x,y) = dhxa(y,x) + dhxb(y,x) + dppija(x,y)
      end do
    end do
    do a = noca + 1, nbf
      do b = noca + 1, a
        dwmo(a,b) = dhxb(b,a)
      end do
    end do

    do i = 1, nbf
      dwmo(i,i) = 0.5_dp*dwmo(i,i)
    end do
    dwmo = -dwmo
  end subroutine differentiate_w_mo_kernel


  function mrsf_w_response_status_message(status) result(message)
    integer, intent(in) :: status
    character(len=96) :: message

    select case (status)
    case (MRSF_W_SUCCESS)
      message = 'success'
    case (MRSF_W_BAD_DIMENSION)
      message = 'inconsistent basis, occupation, packed-response, or perturbation dimensions'
    case (MRSF_W_BAD_REFERENCE)
      message = 'only the two-SOMO ROHF triplet reference is supported'
    case (MRSF_W_BAD_TARGET)
      message = 'only MRSF singlet and triplet target states are supported'
    case (MRSF_W_UMRSF_UNSUPPORTED)
      message = 'UMRSF W response is not verified'
    case (MRSF_W_CAM_UNSUPPORTED)
      message = 'CAM or long-range-corrected W response is not verified'
    case (MRSF_W_ALLOCATION_FAILURE)
      message = 'allocation failed while forming the MRSF W response'
    case default
      message = 'unknown MRSF W response status'
    end select
  end function mrsf_w_response_status_message

end module tdhf_mrsf_hessian_w_response_mod
