!> @file afqmc_trial.F90
!>
!> @brief AFQMC trial wavefunctions built inside OpenQP.
!>
!> The vendored kernels read a multideterminant trial from an `MRSFDAT`/`SFDAT`
!> text file produced by some upstream program.  That still works, but OpenQP
!> computes MRSF-TDDFT itself, so this file turns one of its own MRSF roots into
!> the determinant expansion directly -- no file, no external CI step.
!>
!> Spin-flip determinant basis
!> ---------------------------
!> MRSF works from a high-spin (M_S=+1) ROHF/UHF reference with `nocca` alpha
!> and `noccb` beta occupied orbitals; the two open-shell orbitals are
!> `noccb+1` and `noccb+2 = nocca`.  A response amplitude X(i,j) moves one alpha
!> electron out of occupied orbital `i` into beta orbital `j > noccb`, so it
!> corresponds to exactly one M_S=0 determinant:
!>
!>     alpha occupied : 1..nocca without i          (nocca-1 electrons)
!>     beta  occupied : 1..noccb together with j    (noccb+1 electrons)
!>
!> `mrsfxvec` (tdhf_mrsf_lib) turns the stored, compressed eigenvector into
!> exactly these determinant amplitudes: it is what supplies the spin-complete
!> partner in the open-shell block, i.e. the +/- 1/sqrt(2) pair that
!> distinguishes MRSF from a plain SF vector.  Reusing it is what makes this
!> trial spin complete rather than "an SF vector relabelled as MRSF".

module afqmc_trial_mod

  use precision, only: dp

  implicit none

  private
  public :: build_mrsf_state_trial, read_trial_header

  character(len=*), parameter :: module_name = "afqmc_trial_mod"

contains

!> @brief Build the AFQMC determinant trial from MRSF-TDDFT root @p istate.
!>
!> @param[in]  infos    OpenQP handle (after the MRSF energy step has run)
!> @param[in]  istate   1-based MRSF root
!> @param[in]  thresh   drop determinants whose |coefficient| is below this
!> @param[out] nocc     (alpha, beta) electron counts of the M_S=0 trial
!> @param[out] coeff    determinant coefficients, normalized
!> @param[out] occ      occupied-orbital lists, shape (nbf, 2, ndet)
  subroutine build_mrsf_state_trial(infos, istate, thresh, nocc, coeff, occ)

    use types, only: information
    use io_constants, only: iw
    use messages, only: show_message, WITH_ABORT
    use oqp_tagarray_driver
    use tdhf_mrsf_lib, only: mrsfxvec

    type(information), target, intent(inout) :: infos
    integer, intent(in) :: istate
    real(kind=dp), intent(in) :: thresh
    integer, intent(out) :: nocc(2)
    complex(kind=dp), allocatable, intent(out) :: coeff(:)
    integer, allocatable, intent(out) :: occ(:,:,:)

    real(kind=dp), contiguous, pointer :: bvec_mo(:,:), td_energies(:)
    character(len=*), parameter :: tags_td(2) = (/ character(len=80) :: &
      OQP_td_bvec_mo, OQP_td_energies /)

    real(kind=dp), allocatable :: xexp(:)
    integer, allocatable :: keep(:)
    real(kind=dp) :: amp, norm2, sgn
    integer :: nbf, nocca, noccb, xvec_dim, nstates
    integer :: i, j, ij, idet, ndet, k, ia, ib

    nbf = infos%basis%nbf
    nocca = int(infos%mol_prop%nelec_A)
    noccb = int(infos%mol_prop%nelec_B)

    if (nocca - noccb /= 2) call show_message( &
      "afqmc: an MRSF trial needs the high-spin (M_S=+1) reference, "// &
      "i.e. [scf] multiplicity = 3", WITH_ABORT)

    call data_has_tags(infos%dat, tags_td, module_name, "build_mrsf_state_trial", &
                       WITH_ABORT)
    call tagarray_get_data(infos%dat, OQP_td_bvec_mo, bvec_mo)
    call tagarray_get_data(infos%dat, OQP_td_energies, td_energies)

    xvec_dim = nocca*(nbf - noccb)
    nstates = size(bvec_mo, 2)
    if (size(bvec_mo, 1) /= xvec_dim) call show_message( &
      "afqmc: OQP::td_bvec_mo does not match the MRSF spin-flip dimension", WITH_ABORT)
    if (istate < 1 .or. istate > nstates) call show_message( &
      "afqmc: the requested MRSF state is outside the computed range", WITH_ABORT)

!   Compressed eigenvector -> determinant amplitudes (this is the step that
!   restores the spin-complete open-shell partner).
    allocate(xexp(xvec_dim))
    call mrsfxvec(infos, bvec_mo(:,istate), xexp)

    nocc(1) = nocca - 1
    nocc(2) = noccb + 1

!   Keep only determinants that carry weight; the multideterminant overlap,
!   force bias and local energy all scale linearly in the kept count.
    allocate(keep(xvec_dim))
    ndet = 0
    do ij = 1, xvec_dim
      if (abs(xexp(ij)) > thresh) then
        ndet = ndet + 1
        keep(ndet) = ij
      end if
    end do
    if (ndet == 0) call show_message( &
      "afqmc: the MRSF trial threshold discarded every determinant", WITH_ABORT)

    allocate(coeff(ndet), occ(nbf, 2, ndet))
    occ = 0

    do idet = 1, ndet
      ij = keep(idet)
!     ij = (j - noccb - 1)*nocca + i, with j the beta MO index.
      i = mod(ij - 1, nocca) + 1
      j = (ij - 1)/nocca + noccb + 1
      amp = xexp(ij)

!     Determinant phase. The amplitude multiplies a^dagger_{j beta} a_{i alpha}
!     |Phi0>, while the trial is stored as an ordered column list; writing the
!     alpha columns in ascending order after deleting i costs (-1)**(i-1). The
!     spin-block reordering contributes a factor that is the same for every
!     determinant, so it is an overall phase and is dropped.
      sgn = 1.0_dp
      if (mod(i - 1, 2) == 1) sgn = -1.0_dp
      coeff(idet) = cmplx(sgn*amp, 0.0_dp, kind=dp)

      ia = 0
      do k = 1, nocca
        if (k == i) cycle
        ia = ia + 1
        occ(ia, 1, idet) = k
      end do

      ib = 0
      do k = 1, noccb
        ib = ib + 1
        occ(ib, 2, idet) = k
      end do
      ib = ib + 1
      occ(ib, 2, idet) = j
    end do

    norm2 = sum(abs(coeff)**2)
    coeff = coeff/sqrt(norm2)

    write(iw,'(/1x,"AFQMC trial: MRSF-TDDFT root ",i0," of ",i0)') istate, nstates
    write(iw,'(1x,"MRSF excitation energy  : ",f20.12," Hartree")') td_energies(istate)
    write(iw,'(1x,"Determinants kept       : ",i0," of ",i0," (|c| > ",es9.2,")")') &
      ndet, xvec_dim, thresh
    write(iw,'(1x,"Trial electrons (a/b)   : ",i0,"/",i0)') nocc(1), nocc(2)

    deallocate(xexp, keep)

  end subroutine build_mrsf_state_trial

!> @brief Read only the header of an SFDAT/MRSFDAT trial file.
!>
!> The vendored header reader validates the electron counts against numbers the
!> caller already has; here the file is the authority, because a spin-flip trial
!> lives in a different M_S sector than the SCF reference.
  subroutine read_trial_header(filename, ndet, nalpha, nbeta)
    use messages, only: show_message, WITH_ABORT
    character(len=*), intent(in) :: filename
    integer, intent(out) :: ndet, nalpha, nbeta
    integer :: iu, ios
    open(newunit=iu, file=trim(filename), status='OLD', action='READ', iostat=ios)
    if (ios /= 0) call show_message( &
      "afqmc: cannot open the trial file "//trim(filename), WITH_ABORT)
    read(iu,*,iostat=ios) ndet, nalpha, nbeta
    close(iu)
    if (ios /= 0 .or. ndet < 1 .or. nalpha < 0 .or. nbeta < 0) call show_message( &
      "afqmc: bad header in the trial file "//trim(filename), WITH_ABORT)
  end subroutine read_trial_header

end module afqmc_trial_mod

!> @brief Fill the trial orbital matrices from occupation lists.
!>
!> Split out of AFQMC_Read_MRSFCIS_Trial so an in-memory expansion (see
!> afqmc_trial_mod) reaches AFQMC_Run through exactly the same construction,
!> normalization and dominant-determinant choice as a file-based one.
!> Free-standing rather than in a module because AFQMC_Run is itself a
!> free-standing vendored subroutine.
SUBROUTINE AFQMC_Trial_From_Occupations(NVar,NSpin,NOcc,NDet,Coeff,DetOcc,Phi,IDom)
   USE IO_Files,     ONLY: IW
   USE MPI_Parallel, ONLY: main_rank
   IMPLICIT NONE
   INTEGER,          INTENT(IN)    :: NVar, NSpin, NDet
   INTEGER,          INTENT(IN)    :: NOcc(1:NSpin)
   DOUBLE COMPLEX,   INTENT(INOUT) :: Coeff(1:NDet)
   INTEGER,          INTENT(IN)    :: DetOcc(1:NVar,1:NSpin,1:NDet)
   DOUBLE COMPLEX,   INTENT(OUT)   :: Phi(1:NVar,1:NVar,1:NSpin,1:NDet)
   INTEGER,          INTENT(OUT)   :: IDom
   INTEGER          :: D, S, I, J, Orb
   DOUBLE PRECISION :: Norm2, MaxCoeff

   IF (NSpin /= 2) ERROR STOP "AFQMC trials require explicit alpha and beta blocks"

   Phi = (0.D0,0.D0)
   DO D = 1, NDet
      DO S = 1, NSpin
         DO I = 1, NOcc(S)
            Orb = DetOcc(I,S,D)
            IF (Orb < 1 .OR. Orb > NVar) ERROR STOP "AFQMC trial orbital index out of range"
            DO J = 1, I-1
               IF (DetOcc(J,S,D) == Orb) ERROR STOP "Duplicate orbital in trial determinant"
            ENDDO
            ! Column ordering is part of the determinant phase convention.
            Phi(Orb,I,S,D) = (1.D0,0.D0)
         ENDDO
      ENDDO
   ENDDO

   Norm2 = SUM(ABS(Coeff)**2)
   IF (Norm2 <= TINY(1.D0)) ERROR STOP "All AFQMC trial coefficients are zero"
   Coeff = Coeff/DSQRT(Norm2)

   IDom = 1
   MaxCoeff = ABS(Coeff(1))
   DO D = 2, NDet
      IF (ABS(Coeff(D)) > MaxCoeff) THEN
         IDom = D
         MaxCoeff = ABS(Coeff(D))
      ENDIF
   ENDDO
   IF (main_rank) THEN
      WRITE(IW,'(1X,"AFQMC trial: ",I0," determinants; dominant determinant ",I0)')NDet,IDom
   ENDIF
END SUBROUTINE AFQMC_Trial_From_Occupations

!> @brief Thread-safe wrapper around the walker re-orthogonalisation.
!>
!> AFQMC_Orthogonalize_Walkers_QR_Z needs a copy of the walker orbitals as its
!> input while it overwrites them. The propagation loop kept one shared scratch
!> array for that, which prevents threading the stabilisation step; here the
!> copy is an automatic array, so every thread gets its own.
SUBROUTINE AFQMC_Stabilize_Walker_Z(NVar,NSpin,NOcc,UW,Det,Ovlp)
   IMPLICIT NONE
   INTEGER,        INTENT(IN)    :: NVar, NSpin
   INTEGER,        INTENT(IN)    :: NOcc(1:NSpin)
   DOUBLE COMPLEX, INTENT(INOUT) :: UW(1:NVar,1:NVar,1:NSpin)
   DOUBLE PRECISION, INTENT(OUT) :: Det
   DOUBLE COMPLEX, INTENT(OUT)   :: Ovlp
   DOUBLE COMPLEX                :: Old(1:NVar,1:NVar,1:NSpin)
   Old = UW
   CALL AFQMC_Orthogonalize_Walkers_QR_Z(NVar,NSpin,NOcc,Old,UW,Det,Ovlp)
END SUBROUTINE AFQMC_Stabilize_Walker_Z
