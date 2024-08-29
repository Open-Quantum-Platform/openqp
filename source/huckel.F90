module huckel

  use precision, only: dp
  use oqp_linalg

  implicit none

  private
  public huckel_guess

contains

 subroutine huckel_guess(ovl, orbitals, infos, basis, huckel_basis)

   use constants, only: tol_int
   use types,     only: information
   use messages,  only: show_message, WITH_ABORT
   use mathlib, only: matrix_invsqrt
   use basis_tools, only: basis_set
   use int1, only: basis_overlap
   use guess, only: corresponding_orbital_projection

   implicit none

   type(information), intent(in) :: infos
   type(basis_set), intent(in) :: basis, huckel_basis

   real(kind=dp) :: ovl(*), orbitals(*)
   integer :: nat, i, ok, l0, l0co, nbf, nbf2, nbf_co, nact, ndoc, nproj

   real(kind=dp), allocatable :: scr(:)
   real(kind=dp), allocatable :: q(:)
   real(kind=dp), allocatable :: vec(:,:)
   real(kind=dp), allocatable :: sco(:,:)

   nbf = basis%nbf
   nbf2 = nbf*(nbf+1)/2
   nat = infos%mol_prop%natom

!  Number of orbitals in MINI basis used in Huckel
   nbf_co = huckel_basis%nbf

   allocate(scr(nbf), &
            q(nbf*nbf), &
            vec(nbf_co,nbf_co),     &
            sco(nbf_co,nbf),     &
            stat=ok)
   if (ok/=0) call show_message('Cannot allocate memory', WITH_ABORT)

!  Get overlap between the minimal basis set and the input basis
   call basis_overlap(sco, basis, huckel_basis, tol=log(10.0d0)*tol_int)

   do i = 1, nbf
     sco(:,i) = sco(:,i)*basis%bfnrm(i) * huckel_basis%bfnrm
   end do

!  Determine which orbitals should be projected
   if (infos%control%scftype == 1) then
     ndoc = infos%mol_prop%nelec/2
     nact = 0
   else if (infos%control%scftype >= 2) then
     ndoc = infos%mol_prop%nelec_b
     nact = infos%mol_prop%nelec_a-infos%mol_prop%nelec_b
   end if

!  Extended Huckel calculation in mini basis set
   call huckel_calc(huckel_basis, vec, l0co, nat, infos%atoms%zn, tol_int)

!  Do at most 5 virtuals from the huckel
   nproj = min(l0co,ndoc+nact+5)

!  Get canonical orbitals in input basis space
   call matrix_invsqrt(ovl, q, nbf, qrnk=l0)
   orbitals(1:nbf*nbf) = q(1:nbf*nbf)

!  Project minimal basis set guess onto the input canonical orbitals,
   call corresponding_orbital_projection(vec, sco, orbitals, ndoc, nact, nproj, nbf, nbf_co, l0)

   deallocate(vec, sco)

   call orthogonalize_orbitals(q, ovl, orbitals, nproj, l0, nbf, nbf)

 end subroutine huckel_guess

!> @brief   Extended Huckel calculation in a Huzinaga minimal basis set
 subroutine huckel_calc(basis, vec, l0co, nat, zan, tol_int)
    use eigen, only: diag_symm_full
    use mathlib, only: unpack_matrix
    use messages, only: show_message, WITH_ABORT
    use basis_tools, only: basis_set
    use atomic_structure_m, only: atomic_structure
    use huckel_lut, only: lneg => huckel_lneg
    use guess, only: mksphar
    use int1, only: overlap
    use mathlib, only: orthogonal_transform_sym
!
    implicit none
!
    type(basis_set), intent(in) :: basis
    real(kind=dp) :: vec(:,:), zan(:)
    integer :: l0co, nat, tol_int

    real(kind=dp), parameter :: BITSY=0.05D+00
    real(kind=dp), parameter :: FUDGE = 1.75d+00/2.0d0
!
    real(kind=dp) :: eneg(18)
    integer :: l1co, l2co, l3co, &
               ncore, nval, i, iat, ish, kat, kt, &
               ierr, n, nucz, atype, i0, irow, j0, j, ival
    integer :: ndval(4)
    logical :: notsp


    real(kind=dp), allocatable :: s2(:,:), h2(:,:), wrk2(:,:)
    real(kind=dp), allocatable :: eig(:), h(:), s(:)
    real(kind=dp), allocatable :: TSH(:), Q(:,:), SCR(:,:)
    integer, allocatable :: llim(:), iulim(:)
    logical, allocatable :: core(:)

    l1co = basis%nbf
    l2co = (l1co*l1co+l1co)/2
    l3co = l1co*l1co

    ncore = 0
    nval = 0

    allocate(s2(l1co,l1co), &
             h2(l1co,l1co), &
             eig(l1co), &
             h(l2co), &
             s(l2co), &
             tsh(l3co), &
             q(l1co,l1co), &
             scr(l1co,8), &
             wrk2(l1co,l1co), &
             source=0.0d0)
    allocate(llim(nat), iulim(nat), source=0)
    allocate(core(l1co))

!   set lower and upper basis functions on each atom,
!   counting is done in terms of spherical harmonics.

    iat = 1
    llim(1) = 1
    do ish = 1, basis%nshell
      kat = basis%origin(ish)
      if (kat /= iat) then
        llim(kat)  = iulim(iat)+1
        iulim(kat) = iulim(iat)
        iat = kat
      end if
      kt = basis%am(ish)
      iulim(kat) = iulim(kat) + 2*kt+1
    end do

!   Compute the minimal basis set's overlap matrix
    call overlap(s, basis, log(10.0d0)*tol_int)

!   transform the overlap matrix to spherical harmonic form.
    call mksphar(tsh,l1co,l0co,notsp, basis)

    if (notsp) then
      call orthogonal_transform_sym(l1co, l0co, s, tsh, l1co, h)
      s(1:l2co) = h(1:l2co)
    end if

!   obtain canonical orthonormal MOs -Q- for the minimal basis set.
    call unpack_matrix(s, q)
    call diag_symm_full(1,l0co,q,l1co,eig,ierr)

    do i = 1, l0co
      q(:,i) = q(:,i) / sqrt(eig(i))
    end do

!   construct the extended Huckel operator -H- directly on top of
!   a copy of the overlap -S- in spherical harmonic space.
    call unpack_matrix(s, h2)

    do n = 1, nat
      nucz = int(zan(n))
!     skip dummy atoms
      if (nucz == 0) cycle

      call huckel_get(nucz,eneg,ncore,nval,ndval,atype)

!     set core orbital energies
      i0 = llim(n) - 1
      do i = 1, ncore
        irow = i0+i
        h2(irow,irow) = eneg(lneg(i,atype))
        core(irow) = .true.
      end do

!     set valence orbital energies.
      i0 = llim(n)+ncore-1
      j0 = 0
      do j = 1, nval
        ival = ndval(j)
        do i = 1, ival
          irow = i0+i
          h2(irow,irow) = eneg(lneg(ncore+j0+i,atype))
          core(irow) = .false.
        end do
        i0 = i0+ival
        j0 = j0+ival
      end do

      if (iulim(n)-i0 > 0) then
        call show_message('Huckel: confusion with MINI basis set', WITH_ABORT)
      end if
    end do
!
!   In view of the very large core orbital energies,
!   scale down all core/core and core/valence overlaps,
!   to reduce the amount of mixing of these types.
!   Because of the large range of parameters, one could
!   consider making BITSY a function of the two diagonal
!   elements of H, or different for C/C versus C/V overlaps.
!
    do i = 2, l0co
      do j = 1, i-1
        if (core(i).or.core(j)) h2(j,i) = bitsy*h2(j,i)
      end do
    end do

!   generate the off-diagonal of the extended Huckel operator
    do i = 2, l0co
       do j=1,i-1
          h2(j,i) = FUDGE*h2(j,i)*(h2(i,i)+h2(j,j))
       end do
    end do

    call diag_symm_full(1,l0co,h2,l1co,eig,ierr)
    if (ierr /= 0) call show_message('Huckel MBS diagonalization failure', WITH_ABORT)

!   orthonormalize appropriately.
    call orthogonalize_orbitals(q,s,h2,l0co,l0co,l1co,l1co)

!   backtransform to the Cartesian MBS.
    if (notsp) then
       call dgemm('N', 'N', l1co, l0co, l0co,&
                   1.0d0, tsh, l1co, &
                          h2,  l1co, &
                   0.0d0, vec, l1co)
    else
       vec = h2
    end if

 end subroutine huckel_calc

!> @brief Orthogonalize orbitals
!> @param[in]     q     matrix of 'canonical orbitals', (ndim x l0)
!> @param[in]     s     symmetrix overlap matrix (nbf x nbf), packed
!> @param[in,out] v     orbitals to transform, (ndim x l0)
!> @param[in]     n     defines, how many orbitals from V space to use
!> @param[in]     l0    dimension of the 'canonical orbitals' space
!> @param[in]     nbf    dimension of the AO basis, nbf >= l0 >= n
!> @param[in]     ndim  leading dimension of q and v
!
!> @details Orbital will be computed in three stesp:
!>   1. compute V = Q^T * S * V
!>   2. orthogonalize first `n` vectors from resulting 'V' space
!>   3. back-transform V = Q*V
!
 subroutine orthogonalize_orbitals(q, s, v, n, l0, nbf, ndim)
    use mathlib, only: unpack_matrix
    implicit none

    integer, intent(in) :: n, l0, nbf, ndim
    real(kind=dp), intent(in) :: q(ndim,*), s(*)
    real(kind=dp), intent(inout) :: v(ndim,*)

    real(kind=dp), allocatable :: u(:,:), tmp(:,:), wrk(:)
    real(kind=dp) :: wrksize(1)
    integer :: lwork, info

    call dgeqrf(l0, n, v, ndim, u, wrksize, -1, info)
    lwork = max(int(wrksize(1)), l0)
    allocate(u(nbf,nbf), tmp(nbf,nbf), wrk(lwork))

    ! 1. Compute Q^T * S * V, store in U
    call unpack_matrix(s, u)
    call dsymm('l', 'u', nbf, n, &
               1.0_dp, u, nbf, &
                       v, ndim, &
               0.0_dp, tmp, nbf)
    call dgemm('t', 'n', l0, n, nbf, &
                1.0_dp, q, ndim, &
                        tmp, nbf, &
                0.0_dp, u, ndim)

    ! 2. Orthogonalize orbitals in U
    call dgeqrf(l0, n, u, ndim, tmp, wrk, lwork, info)
    ! The matrix of orthogonal orbitals U is now stored as a product of elementary reflextors

    ! 3. Transform V = Q*U
    v(:,:l0) = q(:,:l0)
    call dormqr('r', 'n', nbf, l0, n, u, ndim, tmp, v, nbf, wrk, lwork, info)

  end subroutine orthogonalize_orbitals

!>    @brief    return Huckel parameters for atom of charge NUCZ
!
!>    @details  parameters are orbital energies, valence shell
!>              info, and sometimes info about how to use a
!>              minimal basis for semicore ECPs.
!
!>    @param[in]  nucz     nuclear charge (atomic number)
!>    @param[out] eneg     list of orbital energies for input atom `nucz`
!>    @param[out] ncore    the number of core _orbitals_
!>    @param[out] nval     the number of valence _shells_
!>    @param[out] ndval    tells how many functions are in each valence shell
!
!>    @note `ncore` and `ndval` count d and f orbitals as containing 6 and 10 functions, respectively.
 subroutine huckel_get(nucz, eneg, ncore, nval, ndval, atype)

    use huckel_lut, only: huckel_eneg, huckel_ncore, huckel_nval, huckel_ndval
    use messages, only: show_message, WITH_ABORT

    implicit none

    integer, intent(in) :: nucz
    real(kind=dp), intent(out) :: eneg(18)
    integer, intent(out) :: ncore, nval, ndval(4), atype


    select case (nucz)
    case (:0)
      ncore=0
      nval=0
    case (1:103)
      eneg  = huckel_eneg(:,nucz)
      ncore = huckel_ncore(nucz)
      nval  = huckel_nval(nucz)
      ndval = huckel_ndval(:,nucz)
    case default
      call show_message("(A,I5)", " Error!  This atom has nuclear charge ", NUCZ)
      call show_message(" Huckel parameters are unavailable past element Lr", WITH_ABORT)
    end select

    select case(nucz)
    case (:57)   ; atype=1
    case (58:71) ; atype=2
    case (72:86) ; atype=3
    case (87:88) ; atype=4
    case (89:90) ; atype=5
    case (91:103); atype=6
    case default
      call show_message("(A,I5)", " Error!  This atom has nuclear charge ", NUCZ)
      call show_message(" Huckel parameters are unavailable past element Lr", WITH_ABORT)
    end select

 end subroutine huckel_get

end module huckel
