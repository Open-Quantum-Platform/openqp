!>  @brief This module contains types and subroutines to manipulate basis set
!>  @details The main goal of this module is to split SP(L) type shells
!>   onto pair of S and P shells. It significantly simplifies code for
!>   one- and two-electron integrals.
!   REVISION HISTORY:
!>  @date -Sep, 2018- Initial release
!>  @author  Vladimir Mironov
module basis_tools
  use iso_fortran_env, only: real64
  use precision, only: dp
  use atomic_structure_m, only: atomic_structure
  use constants, only: BAS_MXANG, ANGULAR_LABEL, NUM_CART_BF, cart_X, cart_Y, cart_Z
  use io_constants, only: IW
  use parallel, only: par_env_t

  implicit none

  type basis_set
    real(real64), dimension(:), allocatable :: &
      ex, & !< Array of primitive Gaussian exponents
      cc, & !< Array of contraction coefficients
      bfnrm !< Array of normalization constants
    integer, dimension(:), allocatable :: &
      g_offset, &  !< Locations of the first Gaussian in shells
      origin, &    !< Tells which atom the shell is centered on
      am, &        !< Array of shell angular momentum
      ncontr, &    !< Array of contraction degrees
      ao_offset, & !< Indices of shells in the total AO basis
      naos,   &         !< Array of shell's AO numbers
      ecp_zn_num   !< number of electrons removed by ecp
    integer :: &
      nshell = 0, &  !< Number of shells in the basis set
      nprim = 0, &   !< Number of primitive Gaussians in the basis set
      nbf = 0, &     !< Number of basis set functions
      mxcontr = 0, & !< Max. contraction degree
      mxam = 0       !< Max. angular momentum among basis set
    type(atomic_structure), pointer :: atoms

    real(real64), allocatable :: at_mx_dist2(:)
    real(real64), allocatable :: prim_mx_dist2(:)
    real(real64), allocatable :: shell_mx_dist2(:)
    real(real64), allocatable :: shell_centers(:, :)


  contains
    procedure, pass(basis) :: from_file
    procedure, pass(basis) :: append
    procedure, pass(basis) :: dump, load
    procedure, pass(basis) :: normalize_primitives
    procedure, pass(basis) :: normalize_contracted
    procedure, pass(basis) :: set_bfnorms
    procedure, pass(basis) :: reserve => omp_sp_reserve
    procedure, private, pass(basis) :: destroy => omp_sp_destroy

    generic :: aoval => compAOv, compAOVg, compAOVgg
    procedure, pass(basis) :: compAOv
    procedure, pass(basis) :: compAOvg
    procedure, pass(basis) :: compAOvgg

    procedure, pass(basis) :: init_shell_centers

    procedure :: set_screening => comp_basis_mxdists

    procedure, pass(basis) :: basis_broadcast
    procedure, pass(basis) :: bf_label

  end type

  private
  public basis_set
  public bas_norm_matrix
  public bas_denorm_matrix

  interface bas_norm_matrix
      module procedure bas_norm_matrix_tr
      module procedure bas_norm_matrix_sq
  end interface
  interface bas_denorm_matrix
      module procedure bas_denorm_matrix_tr
      module procedure bas_denorm_matrix_sq
  end interface

contains

!>  @brief Append a set of shells to the basis set
!
!   PARAMETERS:
!>  @param[in]  atom            basis set for atom, which contains:
!>                 nshells        number of shells to append
!>                 nprims         number of primites to append
!>                 nbfs           number of basis functions to append
!>                 ang(:)         array of shells angular momentums
!>                 ncontract(:)   array of shells contraction degrees
!>                 ex(:)          array of primitive exponent coefficients
!>                 cc(:)          array of primitive contraction coefficients
!>  @param[in]  atom_index      index of added atom in molecule specification
!>  @param[in[  atom_number     number/charge of added atom in molecule specification
!>  @param[in]  newatom         whether to add a new atom or add to the last one
!>  @param[inout]  err         err flag
!
!>  @author  Vladimir Mironov
!
!   REVISION HISTORY:
!>  @date -Sep, 2021- Initial release
  subroutine append(basis, atom, atom_index, atom_number, err, newatom)
    use basis_library, only: atom_basis_t
    use elements, only: ELEMENTS_LONG_NAME
    class(basis_set), intent(inout) :: basis
    class(atom_basis_t), intent(in) :: atom
    integer, intent(in) :: atom_index, atom_number
    logical, optional, intent(in) :: newatom
    logical, intent(inout) :: err

    integer :: i, atom_id, prim_id, bf_id
    logical :: newatom_int

    associate (nshell => basis%nshell &
               , nprim => basis%nprim &
               , nbf => basis%nbf &
               )
!       Check that basis on atom is not emply
      if (atom%nshells == 0) then
        write (iw, '(A,I0,A)') " *** Warning! Element "//trim(ELEMENTS_LONG_NAME(atom_number))//" with index " &
          , atom_index, " does not have basis functions!"
        err = .true.
        return
      end if

      newatom_int = .true.
      if (present(newatom)) newatom_int = newatom

!       Get atom id for new shells
      atom_id = 1
      if (nshell > 0) then
        atom_id = basis%origin(nshell)
        if (newatom_int) atom_id = atom_id+1
      end if

!       Get initial primitive index for new primitives
      prim_id = 1
      if (nshell > 0) prim_id = basis%g_offset(nshell) &
                                +basis%ncontr(nshell)

!       Get initial basis function index for new basis functions
      bf_id = 1
      if (nshell > 0) bf_id = basis%ao_offset(nshell) &
                              +NUM_CART_BF(basis%am(nshell))

!       Set new atom id for new shells
      basis%origin(nshell+1:nshell+atom%nshells) = atom_id

!       Copy contraction degree and angular momentum types
      basis%ncontr(nshell+1:nshell+atom%nshells) = atom%ncontract(:atom%nshells)
      basis%am(nshell+1:nshell+atom%nshells) = atom%ang(:atom%nshells)

!     Update basis max contraction degree
      basis%mxcontr = max(basis%mxcontr, maxval(atom%ncontract(:atom%nshells)))
!     Update basis max angular momentum
      basis%mxam = max(basis%mxam, &
                                             maxval(atom%ang(:atom%nshells)))

!     Compute nbf (used in integral code)
      basis%naos(nshell+1:nshell+atom%nshells) = NUM_CART_BF(atom%ang(:atom%nshells))

      do i = 1, atom%nshells
!           Compute location of the first primitive for every shell
        basis%g_offset(nshell+i) = prim_id
        prim_id = prim_id+atom%ncontract(i)

!           Compute location of the first basis set function for every shell
        basis%ao_offset(nshell+i) = bf_id
        bf_id = bf_id+NUM_CART_BF(atom%ang(i))
      end do

!       Copy primitive exponential and contraction coefficients
      basis%ex(nprim+1:nprim+atom%nprims) = atom%ex(:atom%nprims)
      basis%cc(nprim+1:nprim+atom%nprims) = atom%cc(:atom%nprims)

!       Update counters
      nshell = nshell+atom%nshells
      nprim = nprim+atom%nprims
      nbf = nbf+atom%nbfs

    end associate

  end subroutine

!--------------------------------------------------------------------------------

!>  @brief Allocate arrays in `BASIS_SET` type variable
!
!   PARAMETERS:
!   @param[inout]  basis    basis_set type variable
!
!>  @author  Vladimir Mironov
!
!   REVISION HISTORY:
!>  @date -Sep, 2021- Initial release
  subroutine omp_sp_reserve(basis, num_shell, num_gauss, num_bf)
    class(basis_set), intent(INOUT) :: basis
    integer, intent(in) :: num_shell, num_gauss, num_bf

    if (allocated(basis%ex)) call basis%destroy() ! cleanup

    basis%nshell = 0
    basis%nprim = 0
    basis%nbf = 0

    allocate (basis%ex(num_gauss), source=0.0d0)
    allocate (basis%cc(num_gauss), source=0.0d0)
    allocate (basis%bfnrm(num_bf), source=0.0d0)

    allocate (basis%g_offset(num_shell), source=0)
    allocate (basis%origin(num_shell), source=0)
    allocate (basis%am(num_shell), source=0)
    allocate (basis%ncontr(num_shell), source=0)
    allocate (basis%ao_offset(num_shell), source=0)
    allocate (basis%naos(num_shell), source=0)

  end subroutine

!--------------------------------------------------------------------------------

!>  @brief Dellocate arrays in `BASIS_SET` type variable
!
!   PARAMETERS:
!   @param[inout]  basis    basis_set type variable
!
!>  @author  Vladimir Mironov
!
!   REVISION HISTORY:
!>  @date -Sep, 2021- Initial release
  subroutine omp_sp_destroy(basis)
    class(basis_set), intent(INOUT) :: basis

    deallocate (basis%ex)
    deallocate (basis%cc)
    deallocate (basis%bfnrm)

    deallocate (basis%g_offset)
    deallocate (basis%origin)
    deallocate (basis%am)
    deallocate (basis%ncontr)
    deallocate (basis%ao_offset)
    deallocate (basis%naos)

  end subroutine

!--------------------------------------------------------------------------------

!>  @brief Remove normalization for primitives
!
!>  @author  Vladimir Mironov
!
!   REVISION HISTORY:
!>  @date -Sep, 2021- Initial release
  subroutine normalize_primitives(basis)
    class(basis_set), intent(inout) :: basis
    integer :: ish, ig, ityp, k1, k2
    real(real64) :: ee

    associate (nshell => basis%nshell &
               , g0 => basis%g_offset &
               , am => basis%am &
               , ncontr => basis%ncontr &
               , naos => basis%naos &
               , ex => basis%ex &
               , cc => basis%cc &
               )
      do ish = 1, nshell
        k1 = g0(ish)
        k2 = g0(ish)+ncontr(ish)-1
        ityp = am(ish)

        do ig = k1, k2
          ee = ex(ig)*2.0d0
          cc(ig) = cc(ig)/sqrt(gauss_norm(ee, ityp))
        end do
      end do

    end associate

  end subroutine

!--------------------------------------------------------------------------------

!>  @brief Normalize *contracted* shells
!
!>  @author  Vladimir Mironov
!
!   REVISION HISTORY:
!>  @date -Sep, 2021- Initial release
  subroutine normalize_contracted(basis)
    class(basis_set), intent(inout) :: basis
    integer :: ish, ig, jg, ityp, iat, k1, k2
    real(real64) :: ee, fact, norm
    real(real64), parameter :: THRESH = 1.0d-10
    real(real64), parameter :: WARNAT = 1.0d-6
    character(len=*), parameter :: &
            WRN_ATOM_NORMALIZATION = &
            '(" *** Warning! Atom",I4," shell",I5," type ",A1,&
            &" has normalization",F13.8)'

    associate (nshell => basis%nshell &
               , g0 => basis%g_offset &
               , am => basis%am &
               , ncontr => basis%ncontr &
               , naos => basis%naos &
               , origin => basis%origin &
               , ex => basis%ex &
               , cc => basis%cc &
               )
      do ish = 1, nshell
        k1 = g0(ish)
        k2 = g0(ish)+ncontr(ish)-1
        ityp = am(ish)
        iat = origin(ish)

        fact = 0.0d0
        do ig = k1, k2
          do jg = k1, ig
            ee = ex(ig)+ex(jg)
            norm = cc(ig)*cc(jg)*gauss_norm(ee, ityp)
            if (ig /= jg) norm = 2*norm
            fact = fact+norm
          end do
        end do

        if (fact > THRESH) fact = 1.0d0/sqrt(fact)

        if ((abs(fact-1.0d0) > WARNAT)) then
          write (*, fmt=WRN_ATOM_NORMALIZATION) &
            iat, nshell, ANGULAR_LABEL(ityp+1:ityp+1), fact
        end if

        cc(k1:k2) = cc(k1:k2)*fact
      end do

    end associate

  end subroutine

!--------------------------------------------------------------------------------

!>  @brief Initialize array of basis function normalization factors
!
!   PARAMETERS:
!   @param[inout]  p(:)    array of normalization factors
!
!>  @author  Vladimir Mironov
!
!   REVISION HISTORY:
!>  @date -Sep, 2018- Initial release
  subroutine set_bfnorms(basis)
    use constants, only: shells_pnrm2
    class(basis_set), intent(inout) :: basis

    integer :: mini, maxi, n, i, ang

    do i = 1, basis%nshell
      maxi = basis%naos(i)
      n = basis%ao_offset(i)
      ang = basis%am(i)
      basis%bfnrm(n:n+maxi-1) = shells_pnrm2(1:maxi,ang)
    end do

  end subroutine

!--------------------------------------------------------------------------------

!>  @brief Scale matrix `A` with matrix \f$ P \cdot P^T \f$
!>  @details `A` is a packed square matrix, `P` is a column vector
!
!   PARAMETERS:
!   @param[inout]  a(:)    triangular matrix (dimension `LD`*(`LD`+1)/2)
!   @param[in]     p(:)    vector (dimension `LD`)
!   @param[in]     ld      leading dimension of P and matrix A in unpacked form
!
!>  @author  Vladimir Mironov
!
!   REVISION HISTORY:
!>  @date -Sep, 2018- Initial release
  subroutine bas_norm_matrix_tr(a, p, ld)

    real(real64), intent(INOUT) :: a(:)
    integer, intent(IN) :: ld
    real(real64), allocatable, intent(IN) :: p(:)

    integer :: i, n

!$omp parallel do private(i,n)
    do i = 1, ld
      n = i*(i-1)/2
      a(n+1:n+i) = a(n+1:n+i)*p(i)*p(1:i)
    end do
!$omp end parallel do

  end subroutine

!--------------------------------------------------------------------------------

!>  @brief Scale matrix `A` with matrix \f$ P \cdot P^T \f$
!>  @details `A` is a full square matrix, `P` is a column vector
!
!   PARAMETERS:
!   @param[inout]  a(:,:)  square matrix (dimension `LD`*`LD`)
!   @param[in]     p(:)    vector (dimension `LD`)
!   @param[in]     ld      leading dimension of P and matrix A in unpacked form
!
!>  @author  Vladimir Mironov
!
!   REVISION HISTORY:
!>  @date -Sep, 2018- Initial release
  subroutine bas_norm_matrix_sq(a, p, ld)

    real(real64), intent(INOUT) :: a(:,:)
    integer, intent(IN) :: ld
    real(real64), allocatable, intent(IN) :: p(:)

    integer :: i, n

!$omp parallel do private(i,n)
    do i = 1, ld
      a(:,i) = a(:,i)*p(i)*p(1:i)
    end do
!$omp end parallel do

  end subroutine

!--------------------------------------------------------------------------------

!>  @brief Scale matrix `A` with matrix \f$ 1/P \cdot 1/P^T \f$
!>  @details `A` is a packed square matrix, `P` is a column vector
!
!   PARAMETERS:
!   @param[inout]  a(:)    triangular matrix (dimension `LD`*(`LD`+1)/2)
!   @param[in]     p(:)    vector (dimension `LD`)
!   @param[in]     ld      leading dimension of P and matrix A in unpacked form
!
!>  @author  Vladimir Mironov
!
!   REVISION HISTORY:
!>  @date -Sep, 2018- Initial release
  subroutine bas_denorm_matrix_tr(a, p, ld)

    real(real64), intent(INOUT) :: a(:)
    integer, intent(IN) :: ld
    real(real64), allocatable, intent(INOUT) :: p(:)

    p = 1/p
    call bas_norm_matrix(a, p, ld)
    p = 1/p

  end subroutine

!--------------------------------------------------------------------------------

!>  @brief Scale matrix `A` with matrix \f$ 1/P \cdot 1/P^T \f$
!>  @details `A` is a full square matrix, `P` is a column vector
!
!   PARAMETERS:
!   @param[inout]  a(:,:)  square matrix (dimension `LD`*`LD`)
!   @param[in]     p(:)    vector (dimension `LD`)
!   @param[in]     ld      leading dimension of P and matrix A in unpacked form
!
!>  @author  Vladimir Mironov
!
!   REVISION HISTORY:
!>  @date -Sep, 2018- Initial release
  subroutine bas_denorm_matrix_sq(a, p, ld)

    real(real64), intent(INOUT) :: a(:,:)
    integer, intent(IN) :: ld
    real(real64), allocatable, intent(INOUT) :: p(:)

    p = 1/p
    call bas_norm_matrix(a, p, ld)
    p = 1/p

  end subroutine

  subroutine dump(basis, LU)
    class(basis_set), intent(in) :: basis
    integer, intent(in) :: LU
    integer :: i
    write (LU, *) basis%nshell, basis%nprim, basis%nbf
    do i = 1, basis%nshell
      write (LU, '(*(I8))') basis%g_offset(i) &
        , basis%origin(i) &
        , basis%am(i) &
        , basis%ncontr(i) &
        , basis%ao_offset(i) &
        , basis%naos(i)
    end do
    do i = 1, basis%nprim
      write (LU, '(*(ES23.15))') basis%ex(i), basis%cc(i)
    end do
    do i = 1, basis%nbf
      write (LU, '(A8, ES23.15)') basis%bf_label(i), basis%bfnrm(i)
    end do
  end subroutine dump

  subroutine load(basis, LU)
    class(basis_set), intent(inout) :: basis
    integer, intent(in) :: LU
    integer :: i, nshell, nprim, nbf
    character(8) :: label

    read (LU, *) nshell, nprim, nbf

    call basis%reserve(nshell, nprim, nbf)

    basis%nshell = nshell
    basis%nprim = nprim
    basis%nbf = nbf

    do i = 1, basis%nshell
      read (LU, '(*(I8))') basis%g_offset(i) &
        , basis%origin(i) &
        , basis%am(i) &
        , basis%ncontr(i) &
        , basis%ao_offset(i) &
        , basis%naos(i)
    end do
    do i = 1, basis%nprim
      read (LU, '(*(ES23.15))') basis%ex(i), basis%cc(i)
    end do
    do i = 1, basis%nbf
      read (LU, '(A8, ES23.15)') label, basis%bfnrm(i)
    end do

    basis%mxcontr = maxval(basis%ncontr(1:basis%nshell))
    basis%mxam = maxval(basis%am(1:basis%nshell))
  end subroutine load

  elemental function gauss_norm(e, la) result(res)
    use constants, only: pi
    real(real64), parameter :: &
      NORMS(0:*) = pi*sqrt(pi)*[1.0d+00, 0.5d+00, 0.75d+00, 1.875d+00, &
                         6.5625d+00, 29.53125d+00, 162.421875d+00]
    real(real64), intent(in) :: e
    integer, intent(in) :: la
    real(real64) :: res
    real(real64) :: f
    f = e*sqrt(e)
    res = NORMS(la)/(f*e**la)
  end function

!> @brief Compute AO values in a point
!> @param[in]    basis  atomic basis set
!> @param[in]    ptxyz  coordinates of a point in space
!> @param[out]   naos   number of significant AOs
!> @param[out]   aov    AO values
!> @author Vladimir Mironov
  subroutine compAOv(basis, ptxyz, naos, &
                     aov)

    use precision, only: fp
    implicit none

    class(basis_set) :: basis

    real(kind=fp), intent(in) :: ptxyz(3)
    integer, intent(OUT) :: naos
    real(KIND=fp), contiguous, intent(OUT) :: aov(:)

    real(KIND=fp) :: &
      vexp, vexp1, vexp2, dum
    integer :: &
      k2, loci, &
      ishell, imomfct, ix, iy, iz, ityp
    real(kind=fp) :: dr1(0:10,3), rsqrd
    integer :: i

    dr1(0,:) = 1

    naos = 0
    do ishell = 1, basis%nshell
      associate ( &
        iatm => basis%origin(ishell), &
        offset => basis%ao_offset(ishell), &
        k1 => basis%g_offset(ishell), &
        ncontr => basis%ncontr(ishell), &
        am => basis%am(ishell), &
        maxi => basis%naos(ishell))
!       Angular intermediates for the density and gradient

        dr1(1,:) = ptxyz - basis%atoms%xyz(:3,iatm)
        rsqrd = sum( dr1(1,:)**2 )

        if (rsqrd<=basis%shell_mx_dist2(ishell)) then
          vexp1 = 0.0_fp
          vexp2 = 0.0_fp
          k2 = ncontr-1+k1
          do imomfct = k1, k2
            if (rsqrd > basis%prim_mx_dist2(imomfct)) cycle
            dum = basis%ex(imomfct)*rsqrd
            vexp = exp(-dum)*basis%cc(imomfct)
            vexp1 = vexp1+vexp
            vexp2 = vexp2+2*vexp*basis%ex(imomfct)
          end do

        else
          aov(offset:offset+maxi-1) = 0.0_fp
          cycle
        end if
        naos = naos+1

        select case (am)
        case (0)
          aov(offset) = vexp1
        case (1)
          aov(offset  ) = vexp1*dr1(1, 1)
          aov(offset+1) = vexp1*dr1(1, 2)
          aov(offset+2) = vexp1*dr1(1, 3)
        case DEFAULT
          do i = 2, am
            dr1(i,:) = dr1(i-1,:)*dr1(1,:)
          end do
          loci = offset-1
          do ityp = 1, maxi
            ix = cart_x(ityp,am)
            iy = cart_y(ityp,am)
            iz = cart_z(ityp,am)

!           Compute AO value at a grid point
            aov(loci+ityp) = vexp1*dr1(ix, 1) &
                                  *dr1(iy, 2) &
                                  *dr1(iz, 3)
          end do
        end select
      end associate

    end do

  end subroutine

!> @brief Compute AO values and gradient in a point
!> @param[in]    basis  atomic basis set
!> @param[in]    ptxyz  coordinates of a point in space
!> @param[out]   naos   number of significant AOs
!> @param[out]   aov    AO values
!> @param[out]   aogx   AO gradient, X component
!> @param[out]   aogy   AO gradient, Y component
!> @param[out]   aogz   AO gradient, Z component
!> @author Vladimir Mironov
  subroutine compAOvg(basis, ptxyz, naos, &
                      aov, aogx, aogy, aogz)

    use precision, only: fp
    implicit none

    class(basis_set) :: basis

    real(kind=fp), intent(in) :: ptxyz(3)
    integer, intent(OUT) :: naos

    real(KIND=fp), contiguous, intent(OUT) :: aov(:)
    real(KIND=fp), contiguous, intent(OUT) :: aogx(:), aogy(:), aogz(:)

    integer :: &
      ishell, loci, &!, iatm, mini, maxi, loci0, & !ifct
      ityp, ix, iy, iz

    real(KIND=fp) :: &
      vexp1, vexp2, x, y, z, &
      xm, ym, zm, &
      xp, yp, zp

    real(KIND=fp) :: &
      dum, vexp
    integer :: &
      k2, imomfct
    real(kind=fp) :: dr1(-1:10,3), rsqrd
    integer :: i

    dr1(:-1,:) = 0
    dr1(0,:) = 1

    naos = 0

    do ishell = 1, basis%nshell
      associate ( &
        iatm => basis%origin(ishell), &
        maxi => basis%naos(ishell), &
        offset => basis%ao_offset(ishell), &
        k1 => basis%g_offset(ishell), &
        ncontr => basis%ncontr(ishell), &
        am => basis%am(ishell))

        dr1(1,:) = ptxyz - basis%atoms%xyz(:3,iatm)
        rsqrd = sum( dr1(1,:)**2 )

        if (rsqrd <= basis%shell_mx_dist2(ishell)) then
          vexp1 = 0.0_fp
          vexp2 = 0.0_fp
          k2 = ncontr-1+k1
          do imomfct = k1, k2
            if (rsqrd > basis%prim_mx_dist2(imomfct)) cycle
            dum = basis%ex(imomfct)*rsqrd
            vexp = exp(-dum)*basis%cc(imomfct)
            vexp1 = vexp1+vexp
            vexp2 = vexp2+2*vexp*basis%ex(imomfct)
          end do
        else
          aov(offset:offset+maxi-1) = 0.0_fp
          aogx(offset:offset+maxi-1) = 0.0_fp
          aogy(offset:offset+maxi-1) = 0.0_fp
          aogz(offset:offset+maxi-1) = 0.0_fp
          cycle
        end if
        naos = naos+1

        loci = offset-1
        select case (am)
        case (0)
!           Special fast code for S functions
          aov(offset) = vexp1
          aogx(offset) = -vexp2*dr1(1, 1)
          aogy(offset) = -vexp2*dr1(1, 2)
          aogz(offset) = -vexp2*dr1(1, 3)

        case (1)
!           Special fast code for P functions
          x = dr1(1, 1)
          y = dr1(1, 2)
          z = dr1(1, 3)

          aov(offset  ) = vexp1*x
          aov(offset+1) = vexp1*y
          aov(offset+2) = vexp1*z

          aogx(offset  ) = vexp1-vexp2*x*x
          aogx(offset+1) =      -vexp2*x*y
          aogx(offset+2) =      -vexp2*x*z

          aogy(offset  ) =      -vexp2*x*y
          aogy(offset+1) = vexp1-vexp2*y*y
          aogy(offset+2) =      -vexp2*y*z

          aogz(offset  ) =      -vexp2*x*z
          aogz(offset+1) =      -vexp2*y*z
          aogz(offset+2) = vexp1-vexp2*z*z

        case DEFAULT
          do i = 2, am+1
              dr1(i,:) = dr1(i-1,:)*dr1(1,:)
          end do
          do ityp = 1, maxi
            ix = cart_x(ityp,am)
            iy = cart_y(ityp,am)
            iz = cart_z(ityp,am)

            aov(loci+ityp) = vexp1*dr1(ix, 1) &
                                  *dr1(iy, 2) &
                                  *dr1(iz, 3)

!               Compute gradient AO value at a grid point
!               Gradient is by the electron (not nuclear) coordinates
            x = dr1(ix, 1)
            y = dr1(iy, 2)
            z = dr1(iz, 3)

!               Gradient minus one component
            xm = ix*dr1(ix-1, 1)
            ym = iy*dr1(iy-1, 2)
            zm = iz*dr1(iz-1, 3)

!               Gradient plus one component
            xp = dr1(ix+1, 1)*vexp2
            yp = dr1(iy+1, 2)*vexp2
            zp = dr1(iz+1, 3)*vexp2

            aogx(loci+ityp) = (-xp+vexp1*xm)*y*z
            aogy(loci+ityp) = (-yp+vexp1*ym)*x*z
            aogz(loci+ityp) = (-zp+vexp1*zm)*x*y
          end do
        end select

      end associate

    end do

  end subroutine

!> @brief Compute AO values, and their 1st and 2nd derivatives in a point
!> @param[in]    basis    atomic basis set
!> @param[in]    ptxyz    coordinates of a point in space
!> @param[out]   naos     number of significant AOs
!> @param[out]   aov      AO values
!> @param[out]   aogx     AO gradient, X component
!> @param[out]   aogy     AO gradient, Y component
!> @param[out]   aogz     AO gradient, Z component
!> @param[out]   aog2xx   AO 2nd derivative, XX component
!> @param[out]   aog2yy   AO 2nd derivative, YY component
!> @param[out]   aog2zz   AO 2nd derivative, ZZ component
!> @param[out]   aog2xy   AO 2nd derivative, XY component
!> @param[out]   aog2yz   AO 2nd derivative, YZ component
!> @param[out]   aog2xz   AO 2nd derivative, XZ component
!> @author Vladimir Mironov
  subroutine compAOvgg(basis, ptxyz, naos, &
                       aov, aogx, aogy, aogz, &
                       aog2xx, aog2yy, aog2zz, aog2xy, aog2yz, aog2xz)

    use precision, only: fp
    implicit none

    class(basis_set) :: basis

    real(kind=fp), intent(in) :: ptxyz(3)
    integer, intent(OUT) :: naos

    real(KIND=fp), contiguous, intent(OUT) :: aov(:)
    real(KIND=fp), contiguous, intent(OUT) :: aogx(:), aogy(:), aogz(:)
    real(KIND=fp), contiguous, intent(OUT) :: &
      aog2xx(:), aog2yy(:), aog2zz(:), &
      aog2xy(:), aog2yz(:), aog2xz(:)

    integer :: &
      ishell, loci, &!, iatm, mini, maxi, loci0, & !ifct
      ityp, ix, iy, iz

    real(KIND=fp) :: &
      vexp1, vexp2, vexp3, &
      xm, ym, zm, &
      xp, yp, zp, &
      x, y, z
    real(KIND=fp) :: &
      dx, dy, dz, dxx, dyy, dzz, dxy, dxz, dyz

    real(KIND=fp) :: &
      xmm, ymm, zmm, &
      xpp, ypp, zpp

    real(KIND=fp) :: &
      dum, vexp, tmp(3)
    integer :: &
      k2, imomfct
    real(kind=fp) :: dr1(-2:10,3), rsqrd
    integer :: i

    dr1(:-1,:) = 0
    dr1(0,:) = 1

    naos = 0

    do ishell = 1, basis%nshell
      associate ( &
        iatm => basis%origin(ishell), &
        maxi => basis%naos(ishell), &
        offset => basis%ao_offset(ishell), &
        k1 => basis%g_offset(ishell), &
        ncontr => basis%ncontr(ishell), &
        am => basis%am(ishell))

        dr1(1,:) = ptxyz - basis%atoms%xyz(:3,iatm)
        rsqrd = sum( dr1(1,:)**2 )

        if (rsqrd <= basis%shell_mx_dist2(ishell)) then
          vexp1 = 0.0_fp
          vexp2 = 0.0_fp
          vexp3 = 0.0_fp
          k2 = ncontr-1+k1
          do imomfct = k1, k2
            if (rsqrd > basis%prim_mx_dist2(imomfct)) cycle
            dum = basis%ex(imomfct)*rsqrd
            vexp = exp(-dum)*basis%cc(imomfct)
            vexp1 = vexp1+vexp
            vexp2 = vexp2+2*vexp*basis%ex(imomfct)
            vexp3 = vexp3+4*vexp*basis%ex(imomfct)*basis%ex(imomfct)
          end do
        else
          aov(offset:offset+maxi-1) = 0.0_fp
          aogx(offset:offset+maxi-1) = 0.0_fp
          aogy(offset:offset+maxi-1) = 0.0_fp
          aogz(offset:offset+maxi-1) = 0.0_fp
          aog2xx(offset:offset+maxi-1) = 0.0_fp
          aog2yy(offset:offset+maxi-1) = 0.0_fp
          aog2zz(offset:offset+maxi-1) = 0.0_fp
          aog2xy(offset:offset+maxi-1) = 0.0_fp
          aog2yz(offset:offset+maxi-1) = 0.0_fp
          aog2xz(offset:offset+maxi-1) = 0.0_fp
          cycle
        end if
        naos = naos+1

        loci = offset-1
        select case (am)
        case (0)
!       Special fast code for S functions
          aov(offset) = vexp1

          x = dr1(1, 1)
          y = dr1(1, 2)
          z = dr1(1, 3)

          aogx(offset) = -vexp2*x
          aogy(offset) = -vexp2*y
          aogz(offset) = -vexp2*z

          aog2xx(offset) = vexp3*x*x-vexp2
          aog2yy(offset) = vexp3*y*y-vexp2
          aog2zz(offset) = vexp3*z*z-vexp2
          aog2xy(offset) = vexp3*x*y
          aog2yz(offset) = vexp3*y*z
          aog2xz(offset) = vexp3*x*z

        case (1)
!       Special fast code for P functions
          x = dr1(1, 1)
          y = dr1(1, 2)
          z = dr1(1, 3)

          aov(offset)   = vexp1*x
          aov(offset+1) = vexp1*y
          aov(offset+2) = vexp1*z

          aogx(offset) = vexp1-vexp2*x*x
          aogx(offset+1) = -vexp2*x*y
          aogx(offset+2) = -vexp2*x*z

          aogy(offset) = -vexp2*x*y
          aogy(offset+1) = vexp1-vexp2*y*y
          aogy(offset+2) = -vexp2*y*z

          aogz(offset) = -vexp2*x*z
          aogz(offset+1) = -vexp2*y*z
          aogz(offset+2) = vexp1-vexp2*z*z

          tmp = vexp3*dr1(1, 1:3)*dr1(1, 1:3)-vexp2

          aog2xx(offset  ) = (tmp(1)-2.0_fp*vexp2)*x
          aog2xx(offset+1) = tmp(1)*y
          aog2xx(offset+2) = tmp(1)*z

          aog2yy(offset  ) = tmp(2)*x
          aog2yy(offset+1) = (tmp(2)-2.0_fp*vexp2)*y
          aog2yy(offset+2) = tmp(2)*z

          aog2zz(offset  ) = tmp(3)*x
          aog2zz(offset+1) = tmp(3)*y
          aog2zz(offset+2) = (tmp(3)-2.0_fp*vexp2)*z

          aog2xy(offset  ) = aog2xx(offset+1)
          aog2xy(offset+1) = aog2yy(offset)
          aog2xy(offset+2) = vexp3*x*y*z

          aog2yz(offset  ) = vexp3*x*y*z
          aog2yz(offset+1) = aog2yy(offset+2)
          aog2yz(offset+2) = aog2zz(offset+1)

          aog2xz(offset  ) = aog2xx(offset+2)
          aog2xz(offset+1) = vexp3*x*y*z
          aog2xz(offset+2) = aog2zz(offset)

        case DEFAULT
          do i = 2, am+2
              dr1(i,:) = dr1(i-1,:)*dr1(1,:)
          end do

          do ityp = 1, maxi
            ix = cart_x(ityp,am)
            iy = cart_y(ityp,am)
            iz = cart_z(ityp,am)

!         0 components:
            x = dr1(ix, 1)
            y = dr1(iy, 2)
            z = dr1(iz, 3)

!         +1 components
            xp = dr1(ix+1, 1)
            yp = dr1(iy+1, 2)
            zp = dr1(iz+1, 3)

!         -1 components
            xm = ix*dr1(ix-1, 1)
            ym = iy*dr1(iy-1, 2)
            zm = iz*dr1(iz-1, 3)

!         +2 components
            xpp = dr1(ix+2, 1)
            ypp = dr1(iy+2, 2)
            zpp = dr1(iz+2, 3)

!         -2 components
            xmm = ix*(ix-1)*dr1(ix-2, 1)
            ymm = iy*(iy-1)*dr1(iy-2, 2)
            zmm = iz*(iz-1)*dr1(iz-2, 3)

!         AO value:
            aov(loci+ityp) = vexp1*x*y*z

!         First AO derivatives:
            dx = -vexp2*xp + vexp1*xm
            dy = -vexp2*yp + vexp1*ym
            dz = -vexp2*zp + vexp1*zm
            aogx(loci+ityp) = dx*y*z
            aogy(loci+ityp) = x*dy*z
            aogz(loci+ityp) = x*y*dz

!         Second AO derivatives:
            dxx = vexp3*xpp - vexp2*(2*ix+1)*x + vexp1*xmm
            dyy = vexp3*ypp - vexp2*(2*iy+1)*y + vexp1*ymm
            dzz = vexp3*zpp - vexp2*(2*iz+1)*z + vexp1*zmm
            aog2xx(loci+ityp) = dxx*y*z
            aog2yy(loci+ityp) = x*dyy*z
            aog2zz(loci+ityp) = x*y*dzz

            dxy = vexp3*xp*yp - vexp2*(xp*ym+xm*yp) + vexp1*xm*ym
            dyz = vexp3*yp*zp - vexp2*(yp*zm+ym*zp) + vexp1*ym*zm
            dxz = vexp3*xp*zp - vexp2*(xp*zm+xm*zp) + vexp1*xm*zm
            aog2xy(loci+ityp) = dxy*z
            aog2yz(loci+ityp) = dyz*x
            aog2xz(loci+ityp) = dxz*y

          end do
        end select

      end associate

    end do

  end subroutine

!-------------------------------------------------------------------------------

!> @brief Compute maximum extent of basis set primitives up to a given
!>  tolerance.
!> @details Find the largest root of the eqn.: r**n * exp(-a*r**2) = tol,
!>  and fill prim_mx_dist2 array with correspoinding r**2 values
!>  For n==0 (S shells) the solution is trivial.
!>  For n>0 (P,D,F... shells) the equivalent equation is used:
!>      ln(q)/2 - a*q/n - ln(tol)/n == 0, where q = r**2
!>  The solution of this equation is:
!>      q = -n/(2*a) * W_{-1} (-(2*a/n)*tol**(2/n))
!>  where W_{k} (x) - k-th branch of Lambert W function
!>  Assuming that a << 0.5*tol**(-2/n), W_{-1} (x) can be approximated:
!>      W_{-1} (-x) =  log(x) - log(-log(x))
!>  The assumption holds for reasonable basis sets.
!>  Next, the approximated result is then refined by making 1-2
!>  Newton-Raphson steps.
!>  The error of this approximation is (much) less than 10**(-4) Bohr**2
!>  for typical cases:
!>  a < 10^6, n = (1 to 3) (P to F shells), tol = 10**(-10)
!>  a < 10^3, n = (4 to 6) (G to I shells), tol = 10**(-10)
!  TODO: move to basis set related source file
!> @param[in] basis     basis set variable
!> @param[in] mlogtol   -ln(tol) value
!> @author Vladimir Mironov
  subroutine comp_basis_mxdists(basis, mLogTol)
    use precision, only: dp
    class(basis_set), intent(INOUT) :: basis
    real(KIND=dp), intent(IN) :: mLogTol

    real(KIND=dp) :: tmpLogs(1:BAS_MXANG), logLogTol
    integer :: ish, i

    if (allocated(basis%at_mx_dist2)) deallocate (basis%at_mx_dist2)
    allocate (basis%at_mx_dist2(ubound(basis%atoms%zn,1)), source=-1.0_dp)

    if (allocated(basis%prim_mx_dist2)) deallocate (basis%prim_mx_dist2)
    allocate (basis%prim_mx_dist2(basis%nPrim))

    if (allocated(basis%shell_mx_dist2)) deallocate (basis%shell_mx_dist2)
    allocate (basis%shell_mx_dist2(basis%nShell))

    logLogTol = log(mLogTol)
    tmpLogs = [(logLogTol+2.0*mLogTol/i, i=1, BAS_MXANG)]

    do ish = 1, basis%nShell
      basis%shell_mx_dist2(ish) = 0.0
      do i = basis%g_offset(ish), basis%g_offset(ish)+basis%ncontr(ish)-1
        associate (n => basis%am(ish), &
                   a => basis%ex(i), &
                   r2 => basis%prim_mx_dist2(i), &
                   r2sh => basis%shell_mx_dist2(ish))
          if (n == 0) then
!                   Explicit solution:
            r2 = mLogTol/a
          else if (n < 5) then
!                   Approximate result:
            r2 = 0.5d0*n/a*(tmpLogs(n)-log(a))
!                   One NR step:
            r2 = r2*(1-2*(0.5d0*n*log(r2)-a*r2+mLogTol)/(n-2*a*r2))
          else
!                   Approximate result:
            r2 = 0.5d0*n/a*(tmplogs(n)-log(a))
!                   Two NR steps:
            r2 = r2*(1-2*(0.5d0*n*log(r2)-a*r2+mLogTol)/(n-2*a*r2))
            r2 = r2*(1-2*(0.5d0*n*log(r2)-a*r2+mLogTol)/(n-2*a*r2))
          end if
          r2sh = max(r2, r2sh)
        end associate
      end do
    end do

    do ish = 1, basis%nShell
        associate ( iat => basis%origin(ish) )
            basis%at_mx_dist2(iat) = &
                max(basis%at_mx_dist2(iat), basis%shell_mx_dist2(ish))
        end associate
    end do
  end subroutine

!> @brief Initialize array of shell centers used in electronic integral code
  subroutine init_shell_centers(basis)
    class(basis_set), intent(inout) :: basis

    if (allocated(basis%shell_centers)) deallocate(basis%shell_centers)
    allocate(basis%shell_centers(basis%nshell,3))

    basis%shell_centers(1:basis%nshell,1) = basis%atoms%xyz(1,basis%origin(1:basis%nshell))
    basis%shell_centers(1:basis%nshell,2) = basis%atoms%xyz(2,basis%origin(1:basis%nshell))
    basis%shell_centers(1:basis%nshell,3) = basis%atoms%xyz(3,basis%origin(1:basis%nshell))
  end subroutine

!> @brief   Apply selected basis set library to the molecule
!> @param[inout]   basis        applied basis set to system
!> @param[in]      basis_file   path to basis set in OQP basis set format
!> @param[in]      atoms        atoms in systems (pointer to this variable will be set)
!> @param[out]      err        err flag
  subroutine from_file(basis, basis_file, atoms, err)
    use constants, only: bf_names
    use elements, only: ELEMENTS_ATOMNAME, MAX_ELEMENT_Z
    use messages, only: show_message
    use atomic_structure_m, only: atomic_structure
    use basis_library, only: basis_library_t

    implicit none
    class(basis_set), intent(inout) :: basis
    character(len=*), intent(in) :: basis_file
    type(atomic_structure), target, intent(in) :: atoms
    logical, intent(out) :: err
    type(basis_library_t) :: basis_lib

    integer :: i, nshell, nbasis, ngauss, n
    character(len=4) :: bfl
    character(len=2) :: label
    integer :: num, zn, iatshort

    integer, allocatable :: atom_numbers(:)
!
    err=.false.

    basis%atoms => atoms

    call basis_lib%from_file(basis_file)

    atom_numbers = int(atoms%zn)

!   Compute requried memory space for basis
    call basis_lib%calc_req_storage(atom_numbers, &
                                      nshell, ngauss, nbasis)

!   Allocate basis
    call basis%reserve(nshell, ngauss, nbasis)

!   Add atoms one by one
    nbasis = 0
    do i = 1, size(atom_numbers)
      associate (atom => basis_lib%atoms(atom_numbers(i)))
        call basis%append(atom, i, atom_numbers(i), err)
!   ...... Remove after cleanup ......
        nbasis = nbasis+atom%nbfs
!   ^^^^^^ Remove after cleanup ^^^^^^
      end associate
    end do

!   Set primitive norms
    call basis%normalize_primitives
    call basis%normalize_contracted

!   Set BF norms
    call basis%set_bfnorms

  end subroutine from_file

  function bf_label(basis, bf) result(label)
    use constants, only: bf_names
    use elements, only: ELEMENTS_ATOMNAME, MAX_ELEMENT_Z
    implicit none
    class(basis_set), intent(in) :: basis
    integer, intent(in) :: bf
    character(len=8) :: label
    character(len=2) :: atname
    character(len=2) :: atshort

    integer :: i, zn

    i = bf_to_shell(basis, bf)

    associate( iat => basis%origin(i) &
             , ao  => basis%ao_offset(i) &
             , am  => basis%am(i) &
             , nbf => basis%naos(i) &
             )

      zn = min(int(basis%atoms%zn(iat) + &
              basis%ecp_zn_num(iat)), MAX_ELEMENT_Z)

      atname = ""
      if (zn > 0) atname = ELEMENTS_ATOMNAME(zn)

      write (unit=atshort, fmt='(i2)') mod(iat, 100)

      label = atname // atshort // bf_names(bf-ao+1,am)
    end associate
  end function

  function bf_to_shell(basis, bf) result(res)
    class(basis_set), intent(in) :: basis
    integer, intent(in) :: bf
    integer :: res

    integer, parameter :: mxloop = 100
    integer :: i, h, l

    res = 0
    l = 1
    h = basis%nshell
    if (bf < 1 .or. bf > basis%nbf) return

    if (basis%ao_offset(h) <= bf) then
        res = h
        return
    end if

    do i = 1, mxloop
      res = (l+h)/2
      if (l == h-1) exit
      if (basis%ao_offset(res) > bf) then
        h = res
      else
        l = res
      end if
    end do
    if (i > mxloop) res = 0
  end function

  subroutine basis_broadcast(basis, comm, usempi)
    use iso_c_binding, only: c_bool
    class(basis_set), intent(inout) :: basis
    type(par_env_t) :: pe
    integer, parameter :: int32 = selected_int_kind(9)
    integer(kind=int32) :: comm
    integer :: length
    integer :: i, j
    logical(c_bool), intent(in) :: usempi
    ! Initialize MPI
    call pe%init(comm, usempi)

    length = 1

    ! Broadcast the scalar integers
    call pe%bcast(basis%nshell, length)
    call pe%bcast(basis%nprim, length)

    call pe%bcast(basis%nbf, length)

    call pe%bcast(basis%mxcontr, length)

    call pe%bcast(basis%mxam, length)
    if ( pe%rank /= 0) then
      ! Allocate arrays based on the received sizes (on all processes)
      if (.not. allocated(basis%ex)) allocate(basis%ex(basis%nprim))
      if (.not. allocated(basis%cc)) allocate(basis%cc(basis%nprim))
      if (.not. allocated(basis%bfnrm)) allocate(basis%bfnrm(basis%nbf))
      if (.not. allocated(basis%g_offset)) allocate(basis%g_offset(basis%nshell))
      if (.not. allocated(basis%origin)) allocate(basis%origin(basis%nshell))
      if (.not. allocated(basis%am)) allocate(basis%am(basis%nshell))
      if (.not. allocated(basis%ncontr)) allocate(basis%ncontr(basis%nshell))
      if (.not. allocated(basis%ao_offset)) allocate(basis%ao_offset(basis%nshell))
      if (.not. allocated(basis%naos)) allocate(basis%naos(basis%nshell))
      if (.not. allocated(basis%at_mx_dist2)) allocate(basis%at_mx_dist2(basis%nbf))
      if (.not. allocated(basis%prim_mx_dist2)) allocate(basis%prim_mx_dist2(basis%nprim))
      if (.not. allocated(basis%shell_mx_dist2)) allocate(basis%shell_mx_dist2(basis%nshell))
      if (.not. allocated(basis%shell_centers)) allocate(basis%shell_centers(basis%nshell, 3))
    endif

    ! Broadcast the arrays
    call pe%bcast(basis%ex, basis%nprim)
    call pe%bcast(basis%cc, basis%nprim)
    call pe%bcast(basis%bfnrm, basis%nbf)

    call pe%bcast(basis%g_offset, basis%nshell)

    call pe%bcast(basis%origin, basis%nshell)

    call pe%bcast(basis%am, basis%nshell)

    call pe%bcast(basis%ncontr, basis%nshell)

    call pe%bcast(basis%ao_offset, basis%nshell)

    call pe%bcast(basis%naos, basis%nshell)

    if (.not. allocated(basis%ecp_zn_num)) allocate(basis%ecp_zn_num(maxval(basis%origin)))
    call pe%bcast(basis%ecp_zn_num, maxval(basis%origin))


  end subroutine basis_broadcast

end module
