module grd2

!###############################################################################

  use precision, only: dp
  use constants, only: tol_int
  use io_constants, only: iw
  use basis_tools, only: basis_set
  use grd2_rys, only: grd2_int_data_t, grd2_rys_compute, grd2_rys_hess_compute
  use constants, only: BAS_MXANG
  use int2_compute, only: int2_compute_data_t, ints_exchange

!###############################################################################

  implicit none

!###############################################################################

  character(len=*), parameter :: module_name = "grd2"

!###############################################################################

  character(1), parameter :: bfchars(0:6) = ['S', 'P', 'D', 'F', 'G', 'H', 'I']

!###############################################################################

  type, abstract :: grd2_compute_data_t
    logical :: attenuated = .false.
    real(kind=dp) :: mu = 1.0d99
    real(kind=dp) :: hfscale = 1.0d0
    real(kind=dp) :: hfscale2 = 1.0d0 ! can be used in Responce calculations
    real(kind=dp) :: coulscale = 1.0d0
    integer :: cur_pass = 1
  contains
    procedure(grd2_compute_data_t_init), deferred, pass :: init
    procedure(grd2_compute_data_t_clean), deferred, pass :: clean
    procedure(grd2_compute_data_t_get_density), deferred, pass :: get_density
  end type

!###############################################################################

  abstract interface

    subroutine grd2_compute_data_t_init(this)
      import
      implicit none
      class(grd2_compute_data_t), target, intent(inout) :: this
    end subroutine

    subroutine grd2_compute_data_t_clean(this)
      import
      implicit none
      class(grd2_compute_data_t), target, intent(inout) :: this
    end subroutine

    subroutine grd2_compute_data_t_get_density(this, basis, id, dab, dabmax)
      import
      implicit none
      class(grd2_compute_data_t), target, intent(inout) :: this
      type(basis_set), intent(in) :: basis
      integer, intent(in) :: id(4)
      real(kind=dp), target, intent(out) :: dab(*)
      real(kind=dp), intent(out) :: dabmax
    end subroutine
  end interface

  private
  public :: grd2_driver
  public :: grd2_hess_driver
  public :: grd2_compute_data_t

!###############################################################################

contains

!> @brief The driver for the two electron gradient
  subroutine grd2_driver(infos, basis, de, gcomp, &
                         cam, alpha, beta, mu)

    use types, only: information
    use basis_tools, only: basis_set

    implicit none

    type(information), target, intent(inout) :: infos
    type(basis_set), intent(in) :: basis
    real(kind=dp), intent(inout) :: de(:,:)
    class(grd2_compute_data_t), intent(inout) :: gcomp
    logical, optional, intent(in) :: cam
    real(kind=dp), optional, intent(in) :: alpha, beta, mu

    real(kind=dp), allocatable :: de_internal(:,:)

    logical :: do_cam = .false.

    do_cam = infos%dft%cam_flag
    if (present(cam)) do_cam = cam

    if (do_cam) then
      allocate(de_internal, mold=de)

      gcomp%cur_pass = 1
    ! Regular Coulomb and exchange
      de_internal = 0
      gcomp%attenuated = .false.
      gcomp%coulscale = 1.0d0
      gcomp%hfscale = infos%dft%cam_alpha
      gcomp%hfscale2 = infos%tddft%cam_alpha
      if (present(alpha)) gcomp%hfscale2 = alpha
      call grd2_driver_gen(infos, basis, de_internal, gcomp)
      de = de + de_internal

      gcomp%cur_pass = 2
    ! Short-range exchange:
      de_internal = 0
      gcomp%attenuated = .true.
      gcomp%coulscale = 0.0d0
      gcomp%hfscale = infos%dft%cam_beta
      gcomp%hfscale2 = infos%tddft%cam_beta
      if (present(beta)) gcomp%hfscale2 = beta
      gcomp%mu = infos%dft%cam_mu
      if (present(mu)) gcomp%mu = mu
      call grd2_driver_gen(infos, basis, de_internal, gcomp)
      de = de + de_internal
    else
      ! Only adopt the DFT hybrid mixing here for actual DFT calculations
      ! (hamilton>=20).  For pure Hartree-Fock the caller already set the
      ! correct hfscale (=1.0); infos%dft%hfscale is not meaningful in that
      ! case (it is left at its -1.0 sentinel) and must not clobber it.
      if (infos%control%hamilton >= 20) then
        gcomp%hfscale = infos%dft%hfscale
        gcomp%hfscale2 = infos%tddft%hfscale
      end if
      call grd2_driver_gen(infos, basis, de, gcomp)
    end if

  end subroutine

!> @brief The driver for the two electron gradient
  subroutine grd2_driver_gen(infos, basis, de, gcomp)

    use util, only: measure_time
    use messages, only: show_message, WITH_ABORT
    use types, only: information
    use mod_dft_molgrid, only: dft_grid_t
    use int2_pairs, only: int2_pair_storage, int2_cutoffs_t
    use int2_compute, only: petite_quartet_weight, load_petite_shell_map
    use oqp_tagarray_driver
    use parallel, only: par_env_t

    implicit none

    character(len=*), parameter :: subroutine_name = "grd2_driver_gen"

    type(information), target, intent(inout) :: infos
    type(basis_set), intent(in) :: basis
    class(grd2_compute_data_t), intent(inout) :: gcomp
    real(kind=dp), intent(inout) :: de(:,:)

    real(dp), dimension(:), allocatable :: dab
    real(dp), allocatable :: schwarz_ints(:,:)

    real(kind=dp) :: emu2

    real(kind=dp) :: cutoff, cutoff2, dabcut
    real(kind=dp) :: dabmax, gmax
    real(kind=dp) :: zbig
    integer :: numint, i, ij, skip1, skip2, mpi_ij
    integer :: iok, j, k, l, kl
    integer :: maxnbf, maxl
    integer :: q4, sym_nops
    integer(8), contiguous, pointer :: sym_map(:)

    real(kind=dp) :: rtol, dtol

    integer :: itmp
    character(len=64) :: sval
    integer :: ln
    real(kind=dp) :: cutval
    logical :: lstats

    type(grd2_int_data_t) :: gdat

    type(int2_pair_storage) :: ppairs
    type(int2_cutoffs_t) :: cutoffs

    ! tagarray
    integer(4) :: status

    type(par_env_t) :: pe
    call pe%init(infos%mpiinfo%comm, infos%mpiinfo%usempi)

    if (gcomp%attenuated) then
      emu2 = gcomp%mu**2
    end if

!    `cutoff` is the Schwarz screening cutoff
!    `dabcut` is the two particle density cutoff
!
!   The gradient is evaluated at the CONVERGED density, so the derivative-ERI
!   screening may be loosened relative to the SCF Fock build without changing
!   the converged gradient beyond a controllable tolerance.  Opt-in, default-OFF
!   (env unset => byte-identical to the historic 1.0d-10):
!     OQP_GRAD_CUTOFF - Schwarz block cutoff (default 1.0d-10).  1.0d-8 is the
!                       size-robust opt-in (max|dG| <= ~1e-6 a.u. through the
!                       36-atom systems tested).  1.0d-7 is more aggressive but
!                       max|dG| GROWS with system size and exceeds 1e-5 past
!                       ~18-atom HF (DFT/MRSF, HF-exchange scale <=0.5, tolerate
!                       it to larger sizes).  Looser still is unsafe: derivative
!                       integrals amplify the dropped contributions.  See
!                       GRAD_SCREENING_NOTES.md for the per-size/method table.
    cutoff = 1.0d-10
    call get_environment_variable("OQP_GRAD_CUTOFF", sval, ln)
    if (ln > 0) then
      read(sval,*,iostat=itmp) cutval
      if (itmp == 0 .and. cutval > 0.0_dp) cutoff = cutval
    end if
    cutoff2 = cutoff/2.0d+00

    zbig = maxval(basis%ex)
    dabcut = 1.0d-11
    if (zbig>1.0d+06) dabcut = dabcut/10
    if (zbig>1.0d+07) dabcut = dabcut/10

    dtol = 10.0d0**(-tol_int)
    rtol = log(10.0_dp)*tol_int

!   Opt-in screening diagnostics (auto-enabled when the lever is active).
    lstats = (cutoff /= 1.0d-10)
    call get_environment_variable("OQP_GRAD_STATS", sval, ln)
    if (ln > 0) lstats = (sval(1:1)=='1' .or. sval(1:1)=='y' .or. sval(1:1)=='Y' &
                          .or. sval(1:1)=='t' .or. sval(1:1)=='T')

    call cutoffs%set(&
            cutoff_integral_value=dabcut,&
            cutoff_exp=rtol, &
            cutoff_prefactor_pq=dtol, &
            cutoff_prefactor_p=dtol)
    call ppairs%alloc(basis, cutoffs)
    call ppairs%compute(basis, cutoffs)

    ! integrals for screening
    allocate(schwarz_ints(basis%nshell, basis%nshell))
    if (gcomp%attenuated) then
      call ints_exchange(basis, schwarz_ints, emu2)
    else
      call ints_exchange(basis, schwarz_ints)
    end if

!   Initialize the integral block counters to zero
    skip1 = 0
    skip2 = 0
    numint = 0

!   Optional symmetry petite list (valid: the gradient is linear in the
!   quartets and the SCF density is totally symmetric).
    call load_petite_shell_map(infos, basis%nshell, sym_map, sym_nops)

!   Check maximum angular momentum
    if (basis%mxam>BAS_MXANG) then
      call show_message('gradient integrals programmed up to '&
              //bfchars(BAS_MXANG-1)//' functions', with_abort)
    end if

!   Calculate the largest shell type
    maxnbf = (basis%mxam+1)*(basis%mxam+2)/2

!   Square dtol for use in grd2_rys_compute
    dtol = dtol*dtol

!$omp parallel &
!$omp   private ( &
!$omp   gdat, dab, i, j, k, l, ij, maxl, kl, gmax, dabmax, iok, mpi_ij, q4) &
!$omp   reduction(+:skip1, skip2, numint, de)

     allocate(dab(maxnbf**4))

    call gdat%init(basis%mxam, 1, dtol, dabcut, iok)

!$omp barrier
    if (infos%mpiinfo%usempi) then
       mpi_ij = 0
    end if

    do i = 1, basis%nshell
      do j = 1, i
        ij = i*(i-1)/2+j
        if (ppairs%ppid(1,ij)==0) cycle
        if (infos%mpiinfo%usempi) then
           mpi_ij=mpi_ij+1
           if (mod(mpi_ij, pe%size) /= pe%rank) cycle
        end if

!$omp do schedule(dynamic,4) collapse(2)
        do k = 1, i
          do l = 1, i
          maxl = k
          if (k == i) maxl = j
          if (l > maxl) cycle

            kl = k*(k-1)/2+l
            if (ppairs%ppid(1,kl)==0) cycle

            gmax = schwarz_ints(i,j)*schwarz_ints(k,l)

!           Coarse screening, on just the integral value
            if (gmax<cutoff) then
               skip1 = skip1+1
               cycle
            end if

!           Select centers for derivatives
            call gdat%set_ids(basis,i, j, k, l)

            if (all(gdat%skip(:))) cycle

!           Petite list: keep only the orbit representative; the skeleton
!           gradient is symmetrized (projected) afterwards in pyoqp.
            if (sym_nops > 1) then
              q4 = petite_quartet_weight(sym_map, sym_nops, basis%nshell, i, j, k, l)
              if (q4 == 0) cycle
            else
              q4 = 1
            end if

!           Obtain 2 body density for this shell block
            call gcomp%get_density(basis,gdat%id,dab,dabmax)

!           Fine screening on the weighted contribution (see int2_twoei).
            if (dabmax*gmax*real(q4, dp)<cutoff2) then
               skip2 = skip2+1
               cycle
            end if


!           Evaluate derivative integral, and add to the gradient
            numint = numint+1

            if (gcomp%attenuated) then
              call grd2_rys_compute(gdat, ppairs, dab, dabmax, emu2)
            else
              call grd2_rys_compute(gdat, ppairs, dab, dabmax)
            end if

            de(:,gdat%at) = de(:,gdat%at) + real(q4, dp)*gdat%fd

          end do
        end do
!$omp end do

      end do
    end do

    call gdat%clean()
!$omp end parallel

    call pe%allreduce(skip1, 1)
    call pe%allreduce(skip2, 1)
    call pe%allreduce(numint, 1)
    call pe%allreduce(de, size(de))
!   Finish up the final gradient
!   Project rotational contaminant from gradients
!   call dfinal(1)

    if (lstats) then
      write(iw, fmt="(&
        &/1X,'[grd2] screening: cutoff=',ES9.2,'  pass=',I1,&
        &/1X,'[grd2] blocks coarse/fine skipped ',I12,'/',I12,&
        &'  computed ',I12)") &
        cutoff, gcomp%cur_pass, skip1, skip2, numint
    end if

  end subroutine grd2_driver_gen

!###############################################################################

!> @brief Driver for the analytic two-electron contribution to the Hessian
!>        (skeleton/Hellmann-Feynman part: derivatives of the ERIs contracted
!>        with a fixed two-body density). Mirrors grd2_driver_gen but builds the
!>        per-quartet second-derivative block fd2(3,4,3,4) and scatters it into
!>        the (3*natom,3*natom) Hessian by atom pair.
  recursive subroutine grd2_hess_driver(infos, basis, hess, gcomp, &
                                        cam, alpha, beta, mu)

    use util, only: measure_time
    use messages, only: show_message, WITH_ABORT
    use types, only: information
    use int2_pairs, only: int2_pair_storage, int2_cutoffs_t
    use parallel, only: par_env_t

    implicit none

    type(information), target, intent(inout) :: infos
    type(basis_set), intent(in) :: basis
    class(grd2_compute_data_t), intent(inout) :: gcomp
    real(kind=dp), intent(inout) :: hess(:,:)
    logical, optional, intent(in) :: cam
    real(kind=dp), optional, intent(in) :: alpha, beta, mu

    real(dp), dimension(:), allocatable :: dab
    real(dp), allocatable :: schwarz_ints(:,:)
    real(kind=dp), allocatable :: hess_internal(:,:)

    real(kind=dp) :: emu2
    real(kind=dp) :: cutoff, cutoff2, dabcut
    real(kind=dp) :: dabmax, gmax
    real(kind=dp) :: zbig
    integer :: numint, i, ij, skip1, skip2, mpi_ij
    integer :: iok, j, k, l, kl
    integer :: maxnbf, maxl
    integer :: c1, c2, a1, a2, r0, c0
    real(kind=dp) :: rtol, dtol

    type(grd2_int_data_t) :: gdat
    type(int2_pair_storage) :: ppairs
    type(int2_cutoffs_t) :: cutoffs

    type(par_env_t) :: pe

    logical :: do_cam = .false.

    do_cam = infos%dft%cam_flag
    if (present(cam)) do_cam = cam

    if (do_cam) then
      allocate(hess_internal, mold=hess)

      gcomp%cur_pass = 1
      ! Regular Coulomb and exchange.
      hess_internal = 0.0_dp
      gcomp%attenuated = .false.
      gcomp%coulscale = 1.0_dp
      gcomp%hfscale = infos%dft%cam_alpha
      gcomp%hfscale2 = infos%tddft%cam_alpha
      if (present(alpha)) gcomp%hfscale2 = alpha
      call grd2_hess_driver(infos, basis, hess_internal, gcomp, cam=.false.)
      hess = hess + hess_internal

      gcomp%cur_pass = 2
      ! Short-range exchange.
      hess_internal = 0.0_dp
      gcomp%attenuated = .true.
      gcomp%coulscale = 0.0_dp
      gcomp%hfscale = infos%dft%cam_beta
      gcomp%hfscale2 = infos%tddft%cam_beta
      if (present(beta)) gcomp%hfscale2 = beta
      gcomp%mu = infos%dft%cam_mu
      if (present(mu)) gcomp%mu = mu
      call grd2_hess_driver(infos, basis, hess_internal, gcomp, cam=.false.)
      hess = hess + hess_internal

      deallocate(hess_internal)
      return
    end if

    call pe%init(infos%mpiinfo%comm, infos%mpiinfo%usempi)

    if (infos%control%hamilton >= 20 .and. .not. present(cam)) then
      gcomp%hfscale = infos%dft%hfscale
      gcomp%hfscale2 = infos%tddft%hfscale
    end if

    if (gcomp%attenuated) emu2 = gcomp%mu**2

    cutoff = 1.0d-10
    cutoff2 = cutoff/2.0d+00

    zbig = maxval(basis%ex)
    dabcut = 1.0d-11
    if (zbig>1.0d+06) dabcut = dabcut/10
    if (zbig>1.0d+07) dabcut = dabcut/10

    dtol = 10.0d0**(-tol_int)
    rtol = log(10.0_dp)*tol_int

    call cutoffs%set(&
            cutoff_integral_value=dabcut,&
            cutoff_exp=rtol, &
            cutoff_prefactor_pq=dtol, &
            cutoff_prefactor_p=dtol)
    call ppairs%alloc(basis, cutoffs)
    call ppairs%compute(basis, cutoffs)

    allocate(schwarz_ints(basis%nshell, basis%nshell))
    if (gcomp%attenuated) then
      call ints_exchange(basis, schwarz_ints, emu2)
    else
      call ints_exchange(basis, schwarz_ints)
    end if

    skip1 = 0
    skip2 = 0
    numint = 0

    if (basis%mxam>BAS_MXANG) then
      call show_message('hessian integrals programmed up to '&
              //bfchars(BAS_MXANG-1)//' functions', with_abort)
    end if

    maxnbf = (basis%mxam+1)*(basis%mxam+2)/2
    dtol = dtol*dtol

!$omp parallel &
!$omp   private ( &
!$omp   gdat, dab, i, j, k, l, ij, maxl, kl, gmax, dabmax, iok, mpi_ij, &
!$omp   c1, c2, a1, a2, r0, c0) &
!$omp   reduction(+:skip1, skip2, numint, hess)

    allocate(dab(maxnbf**4))

    call gdat%init(basis%mxam, 2, dtol, dabcut, iok)

!$omp barrier
    if (infos%mpiinfo%usempi) then
       mpi_ij = 0
    end if

    do i = 1, basis%nshell
      do j = 1, i
        ij = i*(i-1)/2+j
        if (ppairs%ppid(1,ij)==0) cycle
        if (infos%mpiinfo%usempi) then
           mpi_ij=mpi_ij+1
           if (mod(mpi_ij, pe%size) /= pe%rank) cycle
        end if

!$omp do schedule(dynamic,4) collapse(2)
        do k = 1, i
          do l = 1, i
          maxl = k
          if (k == i) maxl = j
          if (l > maxl) cycle

            kl = k*(k-1)/2+l
            if (ppairs%ppid(1,kl)==0) cycle

            gmax = schwarz_ints(i,j)*schwarz_ints(k,l)
            if (gmax<cutoff) then
               skip1 = skip1+1
               cycle
            end if

            call gdat%set_ids(basis,i, j, k, l)

            ! All four centers on the same atom -> zero contribution
            if (all(gdat%skip(:))) cycle
            ! Differentiate all four centers explicitly (no TI recovery)
            gdat%skip = .false.

            call gcomp%get_density(basis,gdat%id,dab,dabmax)

            if (dabmax*gmax<cutoff2) then
               skip2 = skip2+1
               cycle
            end if

            numint = numint+1

            if (gcomp%attenuated) then
              call grd2_rys_hess_compute(gdat, ppairs, dab, dabmax, emu2)
            else
              call grd2_rys_hess_compute(gdat, ppairs, dab, dabmax)
            end if

            ! Scatter fd2(a1,c1,a2,c2) into the atom-pair Hessian blocks
            do c1 = 1, 4
              r0 = 3*(gdat%at(c1)-1)
              do c2 = 1, 4
                c0 = 3*(gdat%at(c2)-1)
                do a1 = 1, 3
                  do a2 = 1, 3
                    hess(r0+a1, c0+a2) = hess(r0+a1, c0+a2) &
                                       + gdat%fd2(a1,c1,a2,c2)
                  end do
                end do
              end do
            end do

          end do
        end do
!$omp end do

      end do
    end do

    call gdat%clean()
!$omp end parallel

    call pe%allreduce(skip1, 1)
    call pe%allreduce(skip2, 1)
    call pe%allreduce(numint, 1)
    call pe%allreduce(hess, size(hess))

  end subroutine grd2_hess_driver

!###############################################################################

end module grd2
