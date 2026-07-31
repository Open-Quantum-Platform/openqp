!> @file ccsd_t_energy.F90
!>
!> @brief Driver for the closed-shell CCSD(T) ground-state energy method.
!>
!> CCSD(T) is a post-SCF ground-state correlation method: PyOQP converges the
!> RHF reference first, and this driver adds the coupled-cluster singles and
!> doubles correlation energy plus the perturbative triples correction on top.
!> Dispatched from Python as `[input] method = ccsd(t)`, or as `method = ccsd`
!> to stop after the doubles.
!>
!> The AO integrals come from the shared `int2_compute` engine through a
!> collecting consumer (`cc_ao2mo`), are transformed to the MO basis, and are
!> then consumed by the parallel coupled-cluster kernels in `cc_lib`.
module ccsd_t_energy_mod

  implicit none

  private
  public :: ccsd_t_energy

  character(len=*), parameter :: module_name = "ccsd_t_energy_mod"

contains

  !> C-bound entry point: `[input] method = ccsd(t)` dispatches here.
  subroutine ccsd_t_energy_C(c_handle) bind(C, name="ccsd_t_energy")
    use c_interop, only: oqp_handle_t, oqp_handle_get_info
    use types, only: information
    type(oqp_handle_t) :: c_handle
    type(information), pointer :: inf
    inf => oqp_handle_get_info(c_handle)
    call ccsd_t_energy(inf)
  end subroutine ccsd_t_energy_C

  !> Closed-shell CCSD(T) correlation on the converged RHF reference.
  subroutine ccsd_t_energy(infos)

    use precision, only: dp
    use io_constants, only: iw
    use types, only: information
    use basis_tools, only: basis_set
    use printing, only: print_module_info
    use messages, only: show_message, with_abort
    use int2_compute, only: int2_compute_t
    use cc_ao2mo, only: cc_eri_collect_t, cc_build_mo_blocks
    use cc_lib, only: cc_ccsd_t_energy, cc_options_t
    use parallel, only: par_env_t
    use oqp_tagarray_driver, only: tagarray_get_data, OQP_VEC_MO_A, OQP_E_MO_A

    implicit none

    type(information), target, intent(inout) :: infos

    type(basis_set), pointer :: basis
    type(int2_compute_t) :: int2_driver
    type(cc_eri_collect_t), target :: eri_data
    type(cc_options_t) :: opts
    type(par_env_t) :: pe

    real(kind=dp), contiguous, pointer :: mo_a(:,:), mo_energy_a(:)
    real(kind=dp), allocatable, target :: ao(:,:,:,:)
    real(kind=dp), allocatable :: cmo(:,:), eo(:), ev(:)
    real(kind=dp), allocatable :: oooo(:,:,:,:), ooov(:,:,:,:), oovv(:,:,:,:)
    real(kind=dp), allocatable :: ovov(:,:,:,:), ovvv(:,:,:,:), vvvv(:,:,:,:)

    integer :: nbf, nocc, nfzc, no, nv, nmo, i, p, ok
    real(kind=dp) :: e_ref, e_ccsd, e_t, mem_gb, t_ccsd, t_trip
    logical :: converged, do_t

    open(unit=iw, file=infos%log_filename, position="append")
    call print_module_info('CCSD_T_Energy', 'Computing CCSD(T) ground-state correlation')

    ! ---- reference checks -------------------------------------------------
    if (infos%control%scftype /= 1) then
      close(iw)
      call show_message('CCSD(T) requires a closed-shell RHF reference &
                        &([input] scftype=rhf)', with_abort)
    end if
    if (infos%control%hamilton /= 10) then
      close(iw)
      call show_message('CCSD(T) requires an HF reference; remove the DFT &
                        &functional from [input]', with_abort)
    end if

    basis => infos%basis
    basis%atoms => infos%atoms

    nbf  = basis%nbf
    nocc = int(infos%mol_prop%nelec_a)
    nfzc = int(infos%control%cc_nfzc)

    if (nfzc < 0 .or. nfzc >= nocc) then
      close(iw)
      call show_message('CCSD(T): invalid frozen-core count', with_abort)
    end if

    no  = nocc - nfzc
    nmo = nbf - nfzc
    nv  = nmo - no

    if (no <= 0 .or. nv <= 0) then
      close(iw)
      call show_message('CCSD(T): no correlated occupied or virtual orbitals', with_abort)
    end if

    ! ---- memory guard ------------------------------------------------------
    ! The dominant allocations are the AO tensor (nbf^4), the full MO tensor
    ! built inside the transformation (nmo^4) and the ladder integrals (nv^4).
    mem_gb = ( real(nbf,dp)**4 + 2.0_dp*real(nmo,dp)**4 + real(nv,dp)**4 ) &
             * 8.0_dp / 1.073741824e9_dp
    write(iw,'(/2X,A,I0,A,I0,A,I0)') 'CCSD(T): nbf = ', nbf, &
        ', correlated occ = ', no, ', virt = ', nv
    if (nfzc > 0) write(iw,'(2X,A,I0)') 'CCSD(T): frozen core orbitals = ', nfzc
    write(iw,'(2X,A,F8.2,A)') 'CCSD(T): peak integral storage ~', mem_gb, ' GB'

    if (mem_gb > 64.0_dp) then
      close(iw)
      call show_message('CCSD(T): in-core integral storage exceeds 64 GB; &
                        &use a smaller basis or freeze more core orbitals', with_abort)
    end if

    ! ---- reference orbitals ------------------------------------------------
    call tagarray_get_data(infos%dat, OQP_VEC_MO_A, mo_a)
    call tagarray_get_data(infos%dat, OQP_E_MO_A, mo_energy_a)

    allocate(cmo(nbf,nmo), eo(no), ev(nv), stat=ok)
    if (ok /= 0) call show_message('CCSD(T): cannot allocate MO arrays', with_abort)

    ! Drop the frozen core columns; keep the rest in energy order.
    do p = 1, nmo
      cmo(:,p) = mo_a(:, nfzc+p)
    end do
    do i = 1, no
      eo(i) = mo_energy_a(nfzc+i)
    end do
    do i = 1, nv
      ev(i) = mo_energy_a(nocc+i)
    end do

    ! ---- AO integrals ------------------------------------------------------
    allocate(ao(nbf,nbf,nbf,nbf), stat=ok)
    if (ok /= 0) call show_message('CCSD(T): cannot allocate the AO integral &
                                   &tensor -- system too large for in-core CC', with_abort)
    ao = 0.0_dp

    call int2_driver%init(basis, infos)
    call int2_driver%set_screening()

    eri_data%ao => ao
    eri_data%nbf = nbf
    call int2_driver%run(eri_data)
    call eri_data%clean()
    call int2_driver%clean()

    ! ---- transform and slice ----------------------------------------------
    allocate(oooo(no,no,no,no), ooov(no,no,no,nv), oovv(no,no,nv,nv), &
             ovov(no,nv,no,nv), ovvv(no,nv,nv,nv), vvvv(nv,nv,nv,nv), stat=ok)
    if (ok /= 0) call show_message('CCSD(T): cannot allocate MO integral blocks', with_abort)

    call cc_build_mo_blocks(nbf, nmo, no, cmo, ao, &
                            oooo, ooov, oovv, ovov, ovvv, vvvv)
    deallocate(ao)

    ! ---- coupled cluster ---------------------------------------------------
    do_t = infos%control%cc_triples /= 0
    opts%maxit      = int(infos%control%cc_maxit)
    opts%conv       = infos%control%cc_conv
    opts%ndiis      = int(infos%control%cc_ndiis)
    opts%do_triples = do_t
    opts%verbose    = 1
    opts%iw         = iw

    e_ref = infos%mol_energy%energy

    call pe%init(infos%mpiinfo%comm, infos%mpiinfo%usempi)

    call cc_ccsd_t_energy(no, nv, eo, ev, oooo, ooov, oovv, ovov, ovvv, vvvv, &
                          pe, opts, e_ccsd, e_t, converged, &
                          time_ccsd=t_ccsd, time_triples=t_trip)

    deallocate(oooo, ooov, oovv, ovov, ovvv, vvvv, cmo, eo, ev)

    ! ---- report ------------------------------------------------------------
    write(iw,'(/,2X,60("="))')
    if (do_t) then
      write(iw,'(2X,A)') 'CCSD(T)  (coupled cluster, ground state)'
    else
      write(iw,'(2X,A)') 'CCSD  (coupled cluster, ground state)'
    end if
    write(iw,'(2X,60("="))')
    write(iw,'(2X,A,F20.10)') 'E(reference, SCF)      = ', e_ref
    write(iw,'(2X,A,F20.10)') 'E(CCSD, correlation)   = ', e_ccsd
    write(iw,'(2X,A,F20.10)') 'E(CCSD, total)         = ', e_ref + e_ccsd
    if (do_t) then
      write(iw,'(2X,A,F20.10)') 'E((T), correction)     = ', e_t
      write(iw,'(2X,A,F20.10)') 'E(CCSD(T), correlation)= ', e_ccsd + e_t
      write(iw,'(2X,A,F20.10)') 'E(CCSD(T), total)      = ', e_ref + e_ccsd + e_t
    end if
    write(iw,'(2X,A)') repeat('-', 60)
    write(iw,'(2X,A,F14.2,A)') 'CCSD iterations        = ', t_ccsd, ' s'
    if (do_t) write(iw,'(2X,A,F14.2,A)') '(T) correction         = ', t_trip, ' s'
    write(iw,'(2X,60("="),/)')

    if (.not. converged) then
      close(iw)
      call show_message('CCSD did not converge; increase [cc] maxit', with_abort)
    end if

    infos%mol_energy%energy = e_ref + e_ccsd + e_t
    infos%mol_energy%etot   = e_ref + e_ccsd + e_t

    close(iw)

  end subroutine ccsd_t_energy

end module ccsd_t_energy_mod
