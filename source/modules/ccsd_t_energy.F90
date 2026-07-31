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
    use cc_ao2mo, only: cc_eri_collect_t, cc_build_mo_blocks, cc_packed_length, &
                        cc_build_full_mo
    use cc_uhf_lib, only: cc_uhf_spinorb_build, cc_uhf_spinorb_gb, cc_uhf_ccsd_t
    use mp2_lib, only: semicanonicalize
    use cc_lib, only: cc_ccsd_t_energy, cc_options_t
    use parallel, only: par_env_t
!$  use omp_lib, only: omp_get_max_threads
    use oqp_tagarray_driver, only: tagarray_get_data, OQP_VEC_MO_A, OQP_E_MO_A

    implicit none

    type(information), target, intent(inout) :: infos

    type(basis_set), pointer :: basis
    type(int2_compute_t) :: int2_driver
    type(cc_eri_collect_t), target :: eri_data
    type(cc_options_t) :: opts
    type(par_env_t) :: pe

    real(kind=dp), contiguous, pointer :: mo_a(:,:), mo_energy_a(:)
    real(kind=dp), allocatable, target :: gao(:)
    real(kind=dp), allocatable :: cmo(:,:), eo(:), ev(:)
    real(kind=dp), allocatable :: oooo(:,:,:,:), ooov(:,:,:,:), oovv(:,:,:,:)
    real(kind=dp), allocatable :: ovov(:,:,:,:), ovvv(:,:,:,:), vvvv(:,:,:,:)

    integer :: nbf, nocc, nfzc, no, nv, nmo, i, p, ok
    real(kind=dp) :: e_ref, e_ccsd, e_t, mem_gb, t_ccsd, t_trip
    integer :: nthr, niter
    logical :: open_shell
    logical :: converged, do_t

    open(unit=iw, file=infos%log_filename, position="append")
    call print_module_info('CCSD_T_Energy', 'Computing CCSD(T) ground-state correlation')

    ! ---- reference checks -------------------------------------------------
    open_shell = infos%control%scftype /= 1
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
    ! The dominant allocations are the packed AO integrals (nbf^4/8), the
    ! half-transformed intermediate that lives alongside them (nbf^4/4 at
    ! nmo ~ nbf) and the ladder integrals (nv^4).
    mem_gb = ( real(cc_packed_length(nbf),dp) &
               + 0.25_dp*real(nmo,dp)**2*real(nbf,dp)**2 &
               + real(nv,dp)**4 ) * 8.0_dp / 1.073741824e9_dp
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
    allocate(gao(cc_packed_length(nbf)), stat=ok)
    if (ok /= 0) call show_message('CCSD(T): cannot allocate the AO integral &
                                   &store -- system too large for in-core CC', with_abort)
    gao = 0.0_dp

    call int2_driver%init(basis, infos)
    call int2_driver%set_screening()

    eri_data%g => gao
    eri_data%nbf = nbf
    eri_data%npair = nbf*(nbf+1)/2
    call int2_driver%run(eri_data)
    call eri_data%clean()
    call int2_driver%clean()

    ! ---- transform and slice ----------------------------------------------
    allocate(oooo(no,no,no,no), ooov(no,no,no,nv), oovv(no,no,nv,nv), &
             ovov(no,nv,no,nv), ovvv(no,nv,nv,nv), vvvv(nv,nv,nv,nv), stat=ok)
    if (ok /= 0) call show_message('CCSD(T): cannot allocate MO integral blocks', with_abort)

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

    ! Report the decomposition actually in force.  A silent fallback to one
    ! rank looks exactly like a slow run, so make the parallel width visible.
    nthr = 1
    !$ nthr = omp_get_max_threads()
    write(iw,'(2X,A,I0,A,I0)') 'CCSD(T): MPI ranks = ', int(pe%size), &
                               ', OpenMP threads = ', nthr

    if (open_shell) then
      call run_open_shell()
    else
    call cc_build_mo_blocks(nbf, nmo, no, cmo, gao, &
                            oooo, ooov, oovv, ovov, ovvv, vvvv)
    deallocate(gao)

    call cc_ccsd_t_energy(no, nv, eo, ev, oooo, ooov, oovv, ovov, ovvv, vvvv, &
                          pe, opts, e_ccsd, e_t, converged, &
                          time_ccsd=t_ccsd, time_triples=t_trip)

    deallocate(oooo, ooov, oovv, ovov, ovvv, vvvv)
    end if
    deallocate(cmo, eo, ev)

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

  contains

    !> Open-shell (UHF/ROHF) CCSD(T) via the spin-orbital solver.
    !>
    !> Both spins are semicanonicalised first: a ROHF Fock matrix is not
    !> diagonal in its occ-occ and vir-vir blocks and the coupled-cluster
    !> denominators are undefined until it is.  For UHF it is a no-op.
    subroutine run_open_shell()

      use oqp_tagarray_driver, only: OQP_VEC_MO_B, OQP_FOCK_A, OQP_FOCK_B

      real(kind=dp), contiguous, pointer :: mo_b(:,:), fock_a(:), fock_b(:)
      real(kind=dp), allocatable :: ca_sc(:,:), cb_sc(:,:), ea_sc(:), eb_sc(:)
      real(kind=dp), allocatable :: eri_aa(:,:,:,:), eri_bb(:,:,:,:), eri_ab(:,:,:,:)
      real(kind=dp), allocatable :: gso(:,:,:,:), eso(:), etmp(:)
      integer, allocatable :: ord(:)
      integer :: nso, noa, nob, nocc_so, p, q, i, j, ok2
      real(kind=dp) :: so_gb, tw0, tw1

      noa = int(infos%mol_prop%nelec_a) - nfzc
      nob = int(infos%mol_prop%nelec_b) - nfzc
      nocc_so = noa + nob
      nso = 2*nmo
      so_gb = cc_uhf_spinorb_gb(nmo)

      write(iw,'(2X,A,I0,A,I0,A,I0)') 'CCSD(T): open shell, alpha occ = ', noa, &
          ', beta occ = ', nob, ', spin orbitals = ', nso
      write(iw,'(2X,A,F8.2,A)') 'CCSD(T): spin-orbital storage ~', so_gb, ' GB'

      if (nob < 0 .or. nocc_so <= 0 .or. nso-nocc_so <= 0) then
        close(iw)
        call show_message('CCSD(T): no correlated occupied or virtual spin orbitals', &
                          with_abort)
      end if
      ! The open-shell path stores the full spin-orbital tensor, sixteen times
      ! the spatial one.  Refuse early rather than die in the allocator.
      if (so_gb > 32.0_dp) then
        close(iw)
        call show_message('CCSD(T): open-shell integral storage exceeds 32 GB; &
            &freeze more core orbitals or use a smaller basis', with_abort)
      end if

      ! ROHF is not supported yet.  Semicanonicalisation makes the occ-occ and
      ! vir-vir Fock blocks diagonal, but for ROHF the occupied-virtual block
      ! survives -- and these spin-orbital equations drop exactly those f_ov
      ! terms.  Measured 5e-3 Hartree out on CH2 triplet, which is wrong
      ! quietly; refuse instead until the f_ov terms are carried.
      if (infos%control%scftype == 3) then
        close(iw)
        call show_message('CCSD(T): ROHF reference is not supported yet &
            &(the f_ov terms are missing); use [scf] type=uhf for open-shell &
            &coupled cluster', with_abort)
      end if

      call tagarray_get_data(infos%dat, OQP_VEC_MO_B, mo_b)
      call tagarray_get_data(infos%dat, OQP_FOCK_A, fock_a)
      call tagarray_get_data(infos%dat, OQP_FOCK_B, fock_b)

      allocate(ca_sc(nbf,nbf), cb_sc(nbf,nbf), ea_sc(nbf), eb_sc(nbf), stat=ok2)
      if (ok2 /= 0) call show_message('CCSD(T): open-shell MO alloc failed', with_abort)

      ! Semicanonicalise over the FULL orbital space before dropping frozen
      ! core: the rotation mixes within the occupied block, so it has to see
      ! every orbital it is allowed to mix with.
      call semicanonicalize(nbf, int(infos%mol_prop%nelec_a), mo_a, fock_a, ca_sc, ea_sc)
      call semicanonicalize(nbf, int(infos%mol_prop%nelec_b), mo_b, fock_b, cb_sc, eb_sc)

      allocate(eri_aa(nmo,nmo,nmo,nmo), eri_bb(nmo,nmo,nmo,nmo), &
               eri_ab(nmo,nmo,nmo,nmo), stat=ok2)
      if (ok2 /= 0) call show_message('CCSD(T): open-shell MO integral alloc failed', &
                                      with_abort)
      call cc_build_full_mo(nbf, nmo, ca_sc(:,nfzc+1:nbf), ca_sc(:,nfzc+1:nbf), gao, eri_aa)
      call cc_build_full_mo(nbf, nmo, cb_sc(:,nfzc+1:nbf), cb_sc(:,nfzc+1:nbf), gao, eri_bb)
      call cc_build_full_mo(nbf, nmo, ca_sc(:,nfzc+1:nbf), cb_sc(:,nfzc+1:nbf), gao, eri_ab)
      deallocate(gao)

      allocate(gso(nso,nso,nso,nso), eso(nso), etmp(nso), ord(nso), stat=ok2)
      if (ok2 /= 0) call show_message('CCSD(T): spin-orbital tensor alloc failed', &
                                      with_abort)
      call cc_uhf_spinorb_build(nmo, eri_aa, eri_bb, eri_ab, gso)
      deallocate(eri_aa, eri_bb, eri_ab)

      ! Interleaved spin-orbital energies, then permuted so the occupied ones
      ! come first in ascending energy -- the ordering the solver assumes.
      do p = 1, nmo
        eso(2*p-1) = ea_sc(nfzc+p)
        eso(2*p)   = eb_sc(nfzc+p)
      end do
      do p = 1, nso
        ord(p) = p
      end do
      do i = 1, nso-1
        do j = i+1, nso
          if (eso(ord(j)) < eso(ord(i))) then
            q = ord(i); ord(i) = ord(j); ord(j) = q
          end if
        end do
      end do
      do p = 1, nso
        etmp(p) = eso(ord(p))
      end do
      eso = etmp
      gso = gso(ord,ord,ord,ord)

      call system_clock_seconds(tw0)
      call cc_uhf_ccsd_t(nso, nocc_so, eso, gso, int(infos%control%cc_maxit), &
                         infos%control%cc_conv, do_t, e_ccsd, e_t, converged, niter)
      call system_clock_seconds(tw1)
      t_ccsd = tw1 - tw0
      t_trip = 0.0_dp

      write(iw,'(2X,A,I0)') 'CCSD(T): open-shell iterations = ', niter

      deallocate(gso, eso, etmp, ord, ca_sc, cb_sc, ea_sc, eb_sc)

    end subroutine run_open_shell

    subroutine system_clock_seconds(t)
      real(kind=dp), intent(out) :: t
      integer(8) :: c, r
      call system_clock(c, r)
      t = real(c, kind=dp)/real(r, kind=dp)
    end subroutine system_clock_seconds

  end subroutine ccsd_t_energy

end module ccsd_t_energy_mod
