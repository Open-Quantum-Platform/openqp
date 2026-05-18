! 17 Aug 12 - CHC - Initial program
! This moldule contains initial guess routines.
!
module guess
    use precision,     only: dp
    use oqp_linalg
    implicit none

    private
    public get_ab_initio_density
    public get_ab_initio_orbital
    public corresponding_orbital_projection
    public mksphar

contains

!> @brief  this will calculatte the density matrix
subroutine get_ab_initio_density(alpha_density,alpha_orbital,beta_density,beta_orbital, infos, basis)
   use mathlib, only: orb_to_dens
   use messages, only: show_message, with_abort
   use types, only: information
   use basis_tools, only: basis_set
   implicit none
   type(information), intent(in) :: infos
   type(basis_set), intent(in) :: basis
   integer    :: ok
   integer    :: scftype, nocc, na, nb, nbasis, na2
   real(dp)   :: alpha_density(:), &
                 alpha_orbital(:,:)
   real(dp), optional :: beta_density(:), &
                         beta_orbital(:,:)
   real(dp), dimension(:), allocatable :: occno
!
   scftype = infos%control%scftype
   na   = infos%mol_prop%nelec_a
   nb   = infos%mol_prop%nelec_b
   nbasis = basis%nbf
   nocc = infos%mol_prop%nocc

   allocate(occno(nocc), stat=ok)
   if(ok/=0) call show_message('allocation of occno fails', with_abort)

   occno = 2
   if (scftype/=1) occno = 1

   na2 = nocc
   if (scftype/=1) na2 = na

!  alpha density is always calculated for RHF/ROHF/UHF
   call orb_to_dens(alpha_density,alpha_orbital,occno,na2,nbasis,nbasis)

   if (nb == 0) then
     beta_density = 0
     if (scftype == 2) beta_orbital = 0

   else
     select case (scftype)
     case (1)

     case (2)
       call orb_to_dens(beta_density,beta_orbital,occno,nb,nbasis,nbasis)

     case (3)
       call orb_to_dens(beta_density,alpha_orbital,occno,nb,nbasis,nbasis)
     end select
   end if

end subroutine get_ab_initio_density

!> @brief Solve \f$ F C = \eps S C \f$
 subroutine get_ab_initio_orbital(fock_m, orbitals, orbital_e, q_m)
   use messages,  only: show_message, WITH_ABORT
   use mathlib, only: orthogonal_transform
   use eigen, only: diag_symm_full
   use mathlib, only: unpack_matrix, pack_matrix
   implicit none
!
   real(kind=dp) :: orbital_e(*)
   real(kind=dp) :: fock_m(*)
   real(kind=dp) :: orbitals(:,:), q_m(:,:)

   real(kind=dp), allocatable :: wrk(:,:), f2(:,:), tfock(:,:)
   integer  :: info, ok, nbf

   nbf = ubound(orbitals, 1)
   allocate(wrk(nbf,nbf), f2(nbf,nbf), tfock(nbf, nbf), stat=ok)
   if (ok/=0) call show_message('Cannot allocate memory', WITH_ABORT)

   call unpack_matrix(fock_m, f2)

   ! F = S^{-1/2} * F * S^{-1/2}
   ! can use orthogonal_transform subroutine, because S^{-1/2} is symmetric
   call orthogonal_transform('n', nbf, q_m, f2, tfock, wrk)
   deallocate(f2, wrk)

   ! Get orbitals C' = S^{1/2} C
   call diag_symm_full(1, nbf, tfock, nbf, orbital_e, info)

   ! Recover orbitals in AO basis
   ! C = S^{-1/2} * C'
   call dgemm('n', 'n', nbf, nbf, nbf, &
              1.0_dp, q_m,      nbf, &
                      tfock,    nbf, &
              0.0_dp, orbitals, nbf)

 end subroutine get_ab_initio_orbital

!>
!> @brief Corresponding orbital projection
!> @param[in]     nproj  number of orbitals from `vb` to projected onto
!> @param[in]     l0     number of the first orbitals from the `va` space to project `vb` to
!> @param[in,out] va     `a` orbitals on entry,
!>                       replaced by corresponding orbitals on exit
!> @param[in]     vb     `b` orbitals on entry, unchanged on exit
!> @param[in]     sba    overlap integrals between the two bases,
!
!> @note The corresponding orbitals in the `b` basis are not generated.
!> @note H.F.King, R.E.Stanton, H.Kim, R.E.Wyatt, R.G.Parr J.Chem.Phys. 47, 1936-1941 (1967)

 subroutine corresponding_orbital_projection(vb, sba, va, &
                   ndoc, nact, nproj, nbf, l1co, l0)

  use eigen, only: schmd, diag_symm_full

  implicit none

  real(kind=dp) :: vb(l1co,nproj), sba(l1co,nbf), va(nbf,l0)
  integer :: ndoc, nact, nproj, nbf, l1co, l0

  real(kind=dp), allocatable :: &
        vaco(:,:), rhs(:,:), d(:,:), tmp(:,:), scr(:)
  integer, allocatable :: iwrk(:)
  logical, allocatable :: unused(:)
  integer :: ok

  allocate(vaco(nbf,nbf), &
           rhs(nbf,nbf), &
           d(nbf,l0), &
           tmp(nbf,nbf), &
           scr(nbf), &
           source=0.0_dp, &
           stat=ok)
  allocate(iwrk(nbf), &
           source=0, &
           stat=ok)
  allocate(unused(nbf), &
           source=.false., &
           stat=ok)

! Vb^T * Sba
  call dgemm('t', 'n', nproj, nbf, l1co, &
             1.0_dp, vb, l1co, &
                     sba, l1co, &
             0.0_dp, rhs, nbf)

! 1. Project core space
  call project(va, rhs, 1, ndoc, l0, nbf)
! 2. Project active space
  call project(va, rhs, ndoc+1, ndoc+nact, l0, nbf)
! 3. Project virtual space
  call project(va, rhs, ndoc+nact+1, nproj, l0, nbf)

 contains

  subroutine project(va, rhs, minmo, maxmo, l0, nbf)

    implicit none

    real(kind=dp) :: va(nbf,*), rhs(nbf,*)
    integer :: minmo, maxmo, l0, nbf

    integer :: i, j, ok
    integer :: nrest, nmo

    nrest = l0-minmo+1
    nmo  = maxmo-minmo+1
    if (nrest <= 0) return
    if (nmo == 0) return

!   D = Vb^T * Sab * Va
    call dgemm('n','n', nmo, nrest, nbf, &
               1.0_dp, rhs(minmo,1), nbf, &
                       va(1,minmo),nbf, &
               0.0_dp, d, nbf)

!   Diagonalize D^T * D
!   Eigenvalues run 0 to 1, we want the highest ones first,
!   that is why we use factor -1
!   We don't need eigenvalues
    call dsyrk('u', 't', nrest, nmo, &
               -1.0_dp, d, nbf, &
                0.0_dp, vaco, nbf)
    call diag_symm_full(1, nrest, vaco, nbf, scr, ok)

!   The near zero overlap part of the space gets scrambled
!   in the above diagonalization, we can get a symmetry
!   adapted space by Gram-Schmidt instead.
    call schmd(vaco,nmo,nrest,nbf,scr)

!   Rotate va to corresponding orbital set, a'=a*v,
    call dgemm('n', 'n', nbf, nrest, nrest, &
               1.0_dp, va(1,minmo), nbf, &
                       vaco, nbf, &
               0.0_dp, tmp, nbf)

!   Copy projected orbitals back to end of va
    va(:,minmo:l0) = tmp(:,1:nrest)

!   Get overlap between `b` and the `a'` corresponding set,
!   D = Vb^T * Sab * VaCO
    call dgemm('t', 't', nrest, nmo, nbf, &
               1.0_dp, va(1,minmo), nbf, &
                       rhs(minmo,1), nbf, &
               0.0_dp, tmp, nbf)

!   Permute `a'` set to maximum overlap with `b` set
    unused(:nmo) = .true.
    do i = 1, nmo
      j = maxloc(abs(tmp(:nmo,i)), dim=1, mask=unused)
      iwrk(i) = j
      unused(j) = .false.
    end do

    call reorder_columns(va(1,minmo),iwrk,nmo,nbf)

  end subroutine

 end subroutine corresponding_orbital_projection

!> @brief Reorder a set of molecular orbitals
!> @param[in]     iorder  reordering instructions
!> @param[in,out] v       matrix (ldv,n) to reorder
!> @param[in]     ldv     matrix dimension
!> @param[in]     n       matrix dimension
 subroutine reorder_columns(v, iorder, n, ldv)

  use messages, only: show_message, WITH_ABORT

  implicit none

  real(kind=dp), intent(inout) :: v(ldv,*)
  integer, intent(inout) :: iorder(*)
  integer, intent(in) :: n, ldv

  integer :: i, j, k

  do i = 1, n
    if (.not. any(iorder(1:n) == i)) then
      call show_message("(A,I4,A)", "**** Error, element", i, &
            " is missing from reordering instructions", WITH_ABORT)
    end if
  end do

  do i = 1, n
     j = iorder(i)
     call dswap(ldv, v(:,i), 1, v(:,j), 1)
     do k = i+1, n
        if (iorder(k) == i) iorder(k) = j
     end do
  end do

 end subroutine reorder_columns

!
!> @brief   Generate transformation to spherical harmonics basis
!>
!> @details This is used by the MINI basis during the Huckel guess,
!>          to convert any d or f shells to spherical harmonics.
!>
!> @author  MWS: 4/2013 rewrite returns pure s,p,d,f spherical
      subroutine mksphar(w,l1co,l0co,notsp,basis)
      use messages, only: show_message, WITH_ABORT
      use basis_tools, only: basis_set
!
      implicit none
!
!
      type(basis_set), intent(in) :: basis
      logical :: notsp
      real(kind=dp) :: w(l1co,l1co)
      integer :: l1co, l0co
!
      real(kind=dp), parameter :: RT0304 = SQRT( 3.0D+00/ 4.0D+00)
      real(kind=dp), parameter :: RT0920 = SQRT( 9.0D+00/20.0D+00)
      integer :: sph, ncont, n, cart, ndim, l
!
!     Transform to spherical harmonic basis
      w = 0
      sph = 1
      ncont = 0
      notsp = .false.
      do n = 1, basis%nshell
         notsp = notsp .or. (basis%am(n) > 1)
!
!        S,P,L shell is already spherical harmonics.
         if (basis%am(n) <= 1) then
            cart = basis%ao_offset(n)
            ndim = basis%naos(n)
            do l = 1, ndim
              w(cart+l-1, sph+l-1) = 1.0d0
            enddo
            sph = sph + ndim
         end if
!
         if (basis%am(n) == 2) then
            cart = basis%ao_offset(n)
!           true D orbitals, from eg irrep of Oh
            w(cart  ,sph  ) =  rt0304
            w(cart+1,sph  ) = -rt0304
            w(cart  ,sph+1) = -0.5d0
            w(cart+1,sph+1) = -0.5d0
            w(cart+2,sph+1) =  1.0d0
!           true d orbitals, from t2g irrep of oh
            w(cart+3,sph+2) =  1.0d0
            w(cart+4,sph+3) =  1.0d0
            w(cart+5,sph+4) =  1.0d0
            ncont = ncont+1
            sph = sph+5
         end if
!
         if (basis%am(n) == 3) then
            cart = basis%ao_offset(n)
!           T1U (IN OH) COMBINATIONS
            w(cart  ,sph  ) =  1.0d0
            w(cart+5,sph  ) = -rt0920
            w(cart+7,sph  ) = -rt0920
            w(cart+1,sph+1) =  1.0d0
            w(cart+3,sph+1) = -rt0920
            w(cart+8,sph+1) = -rt0920
            w(cart+2,sph+2) =  1.0d0
            w(cart+4,sph+2) = -rt0920
            w(cart+6,sph+2) = -rt0920
!           t2u (in oh) combinations
            w(cart+5,sph+3) =  rt0304
            w(cart+7,sph+3) = -rt0304
            w(cart+3,sph+4) =  rt0304
            w(cart+8,sph+4) = -rt0304
            w(cart+4,sph+5) =  rt0304
            w(cart+6,sph+5) = -rt0304
!           a2u (in oh) combination
            w(cart+9,sph+6) =  1.0d0
            ncont = ncont+3
            sph = sph+7
         end if
!
         if(basis%am(n) > 3) then
            CALL show_message('MKSPHAR: CALLED FOR BASIS WITH G/H/I AOS', WITH_ABORT)
         end if
      end do

      l0co = l1co - ncont

      end subroutine mksphar
end module guess
