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
    use cc_uhf_lib, only: cc_uhf_spinorb_build, cc_uhf_spinorb_gb, cc_uhf_peak_gb, &
                        cc_uhf_ccsd_t
    use mp2_lib, only: semicanonicalize
    use cc_lib, only: cc_ccsd_t_energy, cc_options_t
    use cholesky_direct, only: cholesky_direct_decompose
    use cholesky_eri, only: cholesky_eri_decompose, cholesky_eri_max_vectors, &
                            cholesky_transform_vv, cholesky_transform_block, &
                            cholesky_assemble_mo_blocks
    use parallel, only: par_env_t
    use memory_info, only: oqp_memory_check, oqp_available_memory_gb, &
                           oqp_mem_str, OQP_MEMORY_SAFETY_FRACTION
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
    real(kind=dp) :: mem_ao, mem_mo, mem_solver, rno, rnv
    real(kind=dp), allocatable :: lvec(:,:), bvv(:,:), boo(:,:), bov(:,:)
    integer :: nchol, maxchol
    real(kind=dp) :: chol_tol, chol_err
    logical :: use_chol, chol_trunc, use_direct
    integer :: chol_npass
    real(kind=dp) :: packed_gb
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
    ! Common to both paths: the packed AO integrals (nbf^4/8) and the
    ! half-transformed intermediate that lives alongside them (nbf^4/4 at
    ! nmo ~ nbf).  The open-shell path has its own, larger accounting below,
    ! so it is charged only for what it actually builds here.
    mem_ao = real(cc_packed_length(nbf),dp) &
             + 0.25_dp*real(nmo,dp)**2*real(nbf,dp)**2

    if (open_shell) then
      mem_gb = mem_ao * 8.0_dp / 1.073741824e9_dp
    else
      ! Closed shell has two competing peaks and the larger one decides.
      !
      ! Transformation: the packed AO store and the half-transformed
      ! intermediate are both live while all six MO blocks are being filled --
      ! not just the ladder integrals.  When no and nv are comparable, oovv,
      ! ovov and ovvv together rival vvvv, so counting nv^4 alone understates
      ! the stage badly enough to admit a job that dies in the allocator.
      !
      ! Solution: the AO store is gone by then, but the six blocks stay and
      ! ccsd_iterate adds roughly a dozen o^2v^2-shaped intermediates, four
      ! o v^3 / v^3 o ones, plus a DIIS history of 2*ndiis amplitude pairs.
      rno = real(no,dp)
      rnv = real(nv,dp)
      mem_mo = rno**4 + rno**3*rnv + 2.0_dp*rno**2*rnv**2 + rno*rnv**3
      if (infos%control%cc_cholesky /= 0) then
        ! The ladder integrals are replaced by the vectors and a per-block
        ! assembly buffer.  nchol is not known until the factorisation runs, so
        ! estimate it from the vectors-per-basis-function ratio observed at
        ! this tolerance; the log prints the actual count afterwards.
        mem_mo = mem_mo + real(nbf*(nbf+1)/2,dp)*15.0_dp*real(nbf,dp) &
                        + rnv**2*15.0_dp*real(nbf,dp)
      else
        mem_mo = mem_mo + rnv**4
      end if
      mem_solver = mem_mo &
                 + 14.0_dp*rno**2*rnv**2 &
                 + 2.0_dp*rno*rnv**3 + 2.0_dp*rnv**3*rno &
                 + 2.0_dp*real(max(int(infos%control%cc_ndiis),0),dp) &
                   * (rno*rnv + rno**2*rnv**2)
      mem_gb = max(mem_ao + mem_mo, mem_solver) * 8.0_dp / 1.073741824e9_dp
    end if
    write(iw,'(/2X,A,I0,A,I0,A,I0)') 'CCSD(T): nbf = ', nbf, &
        ', correlated occ = ', no, ', virt = ', nv
    if (nfzc > 0) write(iw,'(2X,A,I0)') 'CCSD(T): frozen core orbitals = ', nfzc
    ! Refuse against what this machine can actually give, not a constant.
    ! A fixed cap is wrong twice over: it lets a laptop start a job that dies
    ! in the allocator, and it refuses a 500 GB node that had the memory.
    if (.not. open_shell) then
      call oqp_memory_check(mem_gb, 'CCSD(T)', &
          'use a smaller basis or freeze more core orbitals', iw)
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
    ! With Cholesky every MO block comes from the vectors, so the packed store
    ! exists only to produce them.  The direct factorisation produces the same
    ! vectors without it, but sweeps the shell-pair list once per pivot block
    ! and is measurably slower wherever both fit -- so it is chosen on memory,
    ! never on speed.
    use_chol = infos%control%cc_cholesky /= 0
    chol_tol = infos%control%cc_cholesky_tol
    packed_gb = real(cc_packed_length(nbf), dp)*8.0_dp/1.073741824e9_dp
    use_direct = .false.
    if (use_chol .and. .not. open_shell) then
      select case (int(infos%control%cc_cholesky_direct))
      case (1); use_direct = .true.
      case (2); use_direct = .false.
      case default
        if (oqp_available_memory_gb() > 0.0_dp) then
          use_direct = packed_gb > &
              OQP_MEMORY_SAFETY_FRACTION*oqp_available_memory_gb()
        end if
      end select
      if (use_direct) then
        write(iw,'(2X,A)') 'CCSD(T): the packed AO store ('// &
            trim(oqp_mem_str(packed_gb))//') would not fit; factorising directly.'
      end if
    end if

    if (.not. use_direct) then
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
    end if

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
    ! Only the closed-shell solver is distributed; the spin-orbital one is
    ! threaded but not, so claiming its rank count would be claiming a
    ! speed-up that is not there -- every rank would repeat the same work.
    nthr = 1
    !$ nthr = omp_get_max_threads()

    if (open_shell) then
      write(iw,'(2X,A,I0)') 'CCSD(T): OpenMP threads = ', nthr
      if (pe%size > 1) then
        write(iw,'(2X,A,I0,A)') 'CCSD(T): the open-shell solver is not MPI-parallel; ' // &
            'solving on rank 0 of ', int(pe%size), ' and broadcasting the result.'
        write(iw,'(2X,A)') 'CCSD(T): give it threads rather than ranks.'
      end if
      ! Every rank running the same solve would multiply both the peak memory
      ! and the CPU time by the rank count for no gain -- the spin-orbital
      ! kernels take no decomposition.  Solve once and hand the answers out.
      if (pe%rank == 0) then
        call run_open_shell()
      else
        e_ccsd = 0.0_dp; e_t = 0.0_dp; t_ccsd = 0.0_dp; t_trip = 0.0_dp
        niter = 0
        converged = .false.
      end if
      if (pe%size > 1) call bcast_open_shell_results()
    else
      write(iw,'(2X,A,I0,A,I0)') 'CCSD(T): MPI ranks = ', int(pe%size), &
                                 ', OpenMP threads = ', nthr

      ! The ladder integrals (nv^4) are the largest array in the run and the
      ! only consumer of the vvvv block.  Factorising the AO integrals lets the
      ! ladder rebuild them per block from an nv^2 x nchol object instead, so
      ! vvvv is never allocated.  It costs flops -- the assembly is nchol/no^2
      ! times the ladder DGEMM -- which is the trade that turns a job that does
      ! not fit into one that does.



      if (use_chol) then
        maxchol = cholesky_eri_max_vectors(nbf, 20)
        allocate(lvec(nbf*(nbf+1)/2, maxchol), stat=ok)
        if (ok /= 0) call show_message('CCSD(T): cannot allocate Cholesky vectors', &
                                       with_abort)
        if (use_direct) then
          call cholesky_direct_decompose(basis, infos, chol_tol, maxchol, lvec, &
                                         nchol, chol_err, chol_trunc, chol_npass)
          write(iw,'(2X,A,I0,A)') 'CCSD(T): direct factorisation used ', &
              chol_npass, ' integral passes'
        else
          call cholesky_eri_decompose(nbf, gao, chol_tol, maxchol, lvec, nchol, &
                                      chol_err, chol_trunc)
        end if
        write(iw,'(2X,A,I0,A,ES9.2,A,F5.1,A)') &
            'CCSD(T): Cholesky vectors = ', nchol, ' (residual ', chol_err, &
            ', ', real(nchol,dp)/real(nbf,dp), ' per basis function)'
        if (chol_trunc) then
          write(iw,'(2X,A)') 'CCSD(T): WARNING: Cholesky truncated at the vector ' // &
              'cap before reaching the tolerance; the correlation energy is ' // &
              'less accurate than [cc] cholesky_tol requests.'
        end if
        allocate(boo(no*no, nchol), bov(no*nv, nchol), bvv(nv*nv, nchol), stat=ok)
        if (ok /= 0) call show_message('CCSD(T): cannot allocate the MO Cholesky &
                                       &blocks', with_abort)
        call cholesky_transform_block(nbf, no, cmo(:, 1:no), no, cmo(:, 1:no), &
                                      lvec, nchol, boo)
        call cholesky_transform_block(nbf, no, cmo(:, 1:no), nv, cmo(:, no+1:nmo), &
                                      lvec, nchol, bov)
        call cholesky_transform_block(nbf, nv, cmo(:, no+1:nmo), nv, &
                                      cmo(:, no+1:nmo), lvec, nchol, bvv)
        deallocate(lvec)
      end if

      allocate(oooo(no,no,no,no), ooov(no,no,no,nv), oovv(no,no,nv,nv), &
               ovov(no,nv,no,nv), ovvv(no,nv,nv,nv), stat=ok)
      if (ok /= 0) call show_message('CCSD(T): cannot allocate MO integral blocks', with_abort)
      if (use_chol) then
        allocate(vvvv(1,1,1,1), stat=ok)
      else
        allocate(vvvv(nv,nv,nv,nv), stat=ok)
      end if
      if (ok /= 0) call show_message('CCSD(T): cannot allocate the ladder integrals', &
                                     with_abort)

      if (use_chol) then
        ! Every block comes from the vectors, so the packed AO store has no
        ! consumer left -- which is the point: it can be released here, and
        ! with the direct factorisation it is never built at all.
        call cholesky_assemble_mo_blocks(no, nv, nchol, boo, bov, bvv, &
                                         oooo, ooov, oovv, ovov, ovvv)
        deallocate(boo, bov)
      else
        call cc_build_mo_blocks(nbf, nmo, no, cmo, gao, &
                                oooo, ooov, oovv, ovov, ovvv, vvvv)
      end if
      if (allocated(gao)) deallocate(gao)

      if (use_chol) then
        call cc_ccsd_t_energy(no, nv, eo, ev, oooo, ooov, oovv, ovov, ovvv, vvvv, &
                              pe, opts, e_ccsd, e_t, converged, &
                              time_ccsd=t_ccsd, time_triples=t_trip, &
                              bvv=bvv, nchol=nchol)
      else
        call cc_ccsd_t_energy(no, nv, eo, ev, oooo, ooov, oovv, ovov, ovvv, vvvv, &
                              pe, opts, e_ccsd, e_t, converged, &
                              time_ccsd=t_ccsd, time_triples=t_trip)
      end if

      deallocate(oooo, ooov, oovv, ovov, ovvv, vvvv)
      if (allocated(bvv)) deallocate(bvv)
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
      use mathlib, only: unpack_matrix

      real(kind=dp), contiguous, pointer :: mo_b(:,:), fock_a(:), fock_b(:)
      real(kind=dp), allocatable :: ca_sc(:,:), cb_sc(:,:), ea_sc(:), eb_sc(:)
      real(kind=dp), allocatable :: eri_aa(:,:,:,:), eri_bb(:,:,:,:), eri_ab(:,:,:,:)
      real(kind=dp), allocatable :: gso(:,:,:,:), eso(:), etmp(:)
      real(kind=dp), allocatable :: fso(:,:), fov(:,:), fao(:,:), fmo_a(:,:), fmo_b(:,:), scr(:,:)
      integer, allocatable :: ord(:)
      integer :: nso, noa, nob, nocc_so, p, q, i, j, ok2
      integer :: nocc_seen, nvir_seen
      real(kind=dp) :: so_gb, peak_gb

      noa = int(infos%mol_prop%nelec_a) - nfzc
      nob = int(infos%mol_prop%nelec_b) - nfzc
      nocc_so = noa + nob
      nso = 2*nmo

      write(iw,'(2X,A,I0,A,I0,A,I0)') 'CCSD(T): open shell, alpha occ = ', noa, &
          ', beta occ = ', nob, ', spin orbitals = ', nso

      if (nob < 0 .or. nocc_so <= 0 .or. nso-nocc_so <= 0) then
        close(iw)
        call show_message('CCSD(T): no correlated occupied or virtual spin orbitals', &
                          with_abort)
      end if

      ! Peak over the whole path, not just the spin-orbital tensor: the solver
      ! intermediates and the DIIS history are the same order of magnitude, and
      ! a guard that ignores them waves jobs through that then die allocating.
      so_gb = cc_uhf_spinorb_gb(nmo)
      peak_gb = cc_uhf_peak_gb(nmo, nocc_so, int(infos%control%cc_ndiis))
      write(iw,'(2X,A,F8.2,A)') 'CCSD(T): spin-orbital tensor  ~', so_gb, ' GB'
      call oqp_memory_check(peak_gb, 'CCSD(T) open shell', &
          'freeze more core orbitals or use a smaller basis', iw)

      call tagarray_get_data(infos%dat, OQP_VEC_MO_B, mo_b)
      call tagarray_get_data(infos%dat, OQP_FOCK_A, fock_a)
      call tagarray_get_data(infos%dat, OQP_FOCK_B, fock_b)

      allocate(ca_sc(nbf,nbf), cb_sc(nbf,nbf), ea_sc(nbf), eb_sc(nbf), stat=ok2)
      if (ok2 /= 0) call show_message('CCSD(T): open-shell MO alloc failed', with_abort)

      ! Semicanonicalise the CORRELATED window only.  The occupied rotation
      ! mixes core with valence, so rotating the full space and then dropping
      ! the first nfzc columns would correlate a different subspace than the
      ! requested frozen core -- and would drop a different spatial orbital
      ! from each spin, since alpha and beta rotate separately.  Freezing
      ! first keeps the correlated space equal to the span of the reference
      ! orbitals nfzc+1..nbf, which is what frozen-core CC means everywhere
      ! else; the rotation within that space is what makes the denominators
      ! well defined and leaves the energy unchanged.
      call semicanonicalize(nbf, int(infos%mol_prop%nelec_a), mo_a, fock_a, ca_sc, ea_sc, &
                            nfzc=nfzc)
      call semicanonicalize(nbf, int(infos%mol_prop%nelec_b), mo_b, fock_b, cb_sc, eb_sc, &
                            nfzc=nfzc)

      allocate(eri_aa(nmo,nmo,nmo,nmo), eri_bb(nmo,nmo,nmo,nmo), &
               eri_ab(nmo,nmo,nmo,nmo), stat=ok2)
      if (ok2 /= 0) call show_message('CCSD(T): open-shell MO integral alloc failed', &
                                      with_abort)
      call cc_build_full_mo(nbf, nmo, ca_sc(:,nfzc+1:nbf), ca_sc(:,nfzc+1:nbf), gao, eri_aa)
      call cc_build_full_mo(nbf, nmo, cb_sc(:,nfzc+1:nbf), cb_sc(:,nfzc+1:nbf), gao, eri_bb)
      call cc_build_full_mo(nbf, nmo, ca_sc(:,nfzc+1:nbf), cb_sc(:,nfzc+1:nbf), gao, eri_ab)
      deallocate(gao)

      allocate(eso(nso), etmp(nso), ord(nso), stat=ok2)
      if (ok2 /= 0) call show_message('CCSD(T): spin-orbital energy alloc failed', &
                                      with_abort)

      ! Interleaved spin-orbital energies, then permuted so the occupied ones
      ! come first in ascending energy -- the ordering the solver assumes.
      do p = 1, nmo
        eso(2*p-1) = ea_sc(nfzc+p)
        eso(2*p)   = eb_sc(nfzc+p)
      end do
      ! Partition by occupation first, sort by energy second.  Ordering all
      ! spin orbitals by energy and then calling the lowest nocc_so of them
      ! occupied assumes the two spins share an Aufbau boundary, and they need
      ! not -- a beta virtual can lie below an alpha occupied.  Where that
      ! happens the permutation silently exchanges an occupied orbital for a
      ! virtual one and the whole correlation treatment is built on a different
      ! determinant than the SCF converged.  The occupations are known here, so
      ! they decide the partition and energy only orders within it.
      nocc_seen = 0
      nvir_seen = 0
      do p = 1, nmo
        if (p <= noa) then
          nocc_seen = nocc_seen + 1
          ord(nocc_seen) = 2*p-1
        else
          nvir_seen = nvir_seen + 1
          ord(nocc_so+nvir_seen) = 2*p-1
        end if
        if (p <= nob) then
          nocc_seen = nocc_seen + 1
          ord(nocc_seen) = 2*p
        else
          nvir_seen = nvir_seen + 1
          ord(nocc_so+nvir_seen) = 2*p
        end if
      end do

      do i = 1, nocc_so-1
        do j = i+1, nocc_so
          if (eso(ord(j)) < eso(ord(i))) then
            q = ord(i); ord(i) = ord(j); ord(j) = q
          end if
        end do
      end do
      do i = nocc_so+1, nso-1
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

      ! Build the tensor already in the sorted order.  Doing it the other way
      ! round -- build interleaved, then `gso = gso(ord,ord,ord,ord)` -- makes
      ! a full second copy of the biggest array in the run, so the true peak
      ! was twice what the guard above reported.
      allocate(gso(nso,nso,nso,nso), stat=ok2)
      if (ok2 /= 0) call show_message('CCSD(T): spin-orbital tensor alloc failed', &
                                      with_abort)
      call cc_uhf_spinorb_build(nmo, eri_aa, eri_bb, eri_ab, gso, ord=ord)
      deallocate(eri_aa, eri_bb, eri_ab)

      ! Occupied-virtual Fock block in the same spin-orbital basis.  It is zero
      ! for UHF but not for ROHF: semicanonicalisation only diagonalises the
      ! occ-occ and vir-vir blocks, and the solver needs what is left over --
      ! both in the amplitude equations and in the correlation energy, where
      ! Brillouin no longer holds.
      allocate(fao(nbf,nbf), fmo_a(nbf,nbf), fmo_b(nbf,nbf), scr(nbf,nbf), &
               fso(nso,nso), stat=ok2)
      if (ok2 /= 0) call show_message('CCSD(T): open-shell Fock alloc failed', with_abort)
      call unpack_matrix(fock_a, fao, 'U')
      call dgemm('n','n', nbf, nbf, nbf, 1.0_dp, fao, nbf, ca_sc, nbf, 0.0_dp, scr, nbf)
      call dgemm('t','n', nbf, nbf, nbf, 1.0_dp, ca_sc, nbf, scr, nbf, 0.0_dp, fmo_a, nbf)
      call unpack_matrix(fock_b, fao, 'U')
      call dgemm('n','n', nbf, nbf, nbf, 1.0_dp, fao, nbf, cb_sc, nbf, 0.0_dp, scr, nbf)
      call dgemm('t','n', nbf, nbf, nbf, 1.0_dp, cb_sc, nbf, scr, nbf, 0.0_dp, fmo_b, nbf)

      fso = 0.0_dp
      do p = 1, nmo
        do q = 1, nmo
          fso(2*p-1, 2*q-1) = fmo_a(nfzc+p, nfzc+q)
          fso(2*p,   2*q  ) = fmo_b(nfzc+p, nfzc+q)
        end do
      end do
      fso = fso(ord,ord)

      allocate(fov(nocc_so, nso-nocc_so), stat=ok2)
      if (ok2 /= 0) call show_message('CCSD(T): fov alloc failed', with_abort)
      fov = fso(1:nocc_so, nocc_so+1:nso)
      deallocate(fao, fmo_a, fmo_b, scr, fso)

      call cc_uhf_ccsd_t(nso, nocc_so, eso, gso, int(infos%control%cc_maxit), &
                         infos%control%cc_conv, do_t, e_ccsd, e_t, converged, niter, &
                         fov=fov, time_ccsd=t_ccsd, time_triples=t_trip, &
                         ndiis_in=int(infos%control%cc_ndiis))

      write(iw,'(2X,A,I0)') 'CCSD(T): open-shell iterations = ', niter

      deallocate(gso, eso, etmp, ord, ca_sc, cb_sc, ea_sc, eb_sc, fov)

    end subroutine run_open_shell

    !> Hand the rank-0 solve's results to every other rank, so all ranks leave
    !> with the same energies and the same view of whether it converged.
    subroutine bcast_open_shell_results()
      real(kind=dp) :: buf(4)
      integer :: iconv(2)

      buf = [e_ccsd, e_t, t_ccsd, t_trip]
      call pe%bcast(buf, 4)
      e_ccsd = buf(1); e_t = buf(2); t_ccsd = buf(3); t_trip = buf(4)

      iconv(1) = merge(1, 0, converged)
      iconv(2) = niter
      call pe%bcast(iconv, 2)
      converged = iconv(1) /= 0
      niter = iconv(2)
    end subroutine bcast_open_shell_results


  end subroutine ccsd_t_energy

end module ccsd_t_energy_mod
