module resp_mod

  use oqp_linalg

  implicit none

  character(len=*), parameter :: module_name = "resp_mod"

  private
  public oqp_resp_charges
  public add_atom_grid
  public resp_charges_C

!--------------------------------------------------------------------------------

contains

!--------------------------------------------------------------------------------

!> @brief Compute ESP charges fitted by modified Merz-Kollman method
!> @details This is a C interface for a Fortran subroutine
!>
!> @param[in]      c_handle     OQP handle
!
!> @author   Vladimir Mironov
!
!     REVISION HISTORY:
!> @date _Mar, 2023_ Initial release
  subroutine resp_charges_C(c_handle) bind(C, name="resp_charges")
    use c_interop, only: oqp_handle_t, oqp_handle_get_info, oqp_handle_refresh_ptr
    use strings, only: Cstring
    use types, only: information
    type(oqp_handle_t) :: c_handle
    type(information), pointer :: inf
    inf => oqp_handle_get_info(c_handle)
    call oqp_resp_charges(inf)
  end subroutine resp_charges_C

!--------------------------------------------------------------------------------

!> @brief Compute ESP charges fitted by modified Merz-Kollman method
!> @note  see refs:
!>        B.Besler, K. Merz, P.Kollman, J. Comput. Chem., 11, 4, 431-439 (1990)
!>        C.Bayly et. al., J. Phys. Chem., 97, 10269-10280 (1993)
!>
!> @param[in]      infos        OQP handle
!
!> @author   Vladimir Mironov
!
!     REVISION HISTORY:
!> @date _Mar, 2023_ Initial release
  subroutine oqp_resp_charges(infos)
    use precision, only: dp
    use io_constants, only: iw
    use oqp_tagarray_driver
    use basis_tools, only: basis_set
    use messages, only: show_message, with_abort
    use types, only: information
    use strings, only: Cstring, fstring
    use mathlib, only: traceprod_sym_packed, triangular_to_full
    use int1, only: electrostatic_potential
    use lebedev, only: lebedev_get_grid
    use elements, only: ELEMENTS_VDW_RADII

    implicit none

    character(len=*), parameter :: subroutine_name = "oqp_resp_charges"

    type(information), target, intent(inout) :: infos

    integer :: nbf, nbf2, ok
    logical :: urohf
    type(basis_set), pointer :: basis
    real(kind=dp), allocatable :: xyz(:,:), wt(:), pot(:)
    real(kind=dp), allocatable :: leb(:,:), lebw(:)
    real(kind=dp), allocatable :: vdwrad(:)
    real(kind=dp), allocatable :: chg(:)
    integer, allocatable :: neigh(:)
    real(kind=dp) :: rms
    real(kind=dp), allocatable :: den(:)
    real(kind=dp), allocatable :: q0(:), alpha(:)

    integer :: nat, npt, nptcur, nadd, nleb
    integer :: i, layer

    integer, parameter :: nlayers = 4
    real(kind=dp), parameter :: &
      layers(nlayers) = [1.4, 1.6, 1.8, 2.0]
    integer, parameter :: npt_layer(nlayers) = [132,152,192,350]
    integer, parameter :: typ_layer(nlayers) = [3,3,3,0]

    logical :: restr

    ! tagarray
    real(kind=dp), contiguous, pointer :: dmat_a(:), dmat_b(:)
    character(len=*), parameter :: tags_alpha(1) = (/ character(len=80) :: &
      OQP_DM_A /)
    character(len=*), parameter :: tags_beta(2) = (/ character(len=80) :: &
      OQP_DM_A, OQP_DM_B /)

    open (unit=IW, file=infos%log_filename, position="append")

    select case(infos%control%esp)
    case(0,1)
      restr = .false.
      write(iw,'(2/)')
      write(iw,'(4x,a)') '======================='
      write(iw,'(4x,a)') 'ESP charges calculation'
      write(iw,'(4x,a)') '======================='
    case (2)
      restr = .true.
      allocate(q0(nat), alpha(nat), source=0.0d0)
      alpha = infos%control%resp_constr
      select case(infos%control%resp_target)
      case(0)
        q0 = 0
      !case(1)
      ! set Mulliken charges
      case default
        error stop 'Unknown RESP charges target'
      end select
      write(iw,'(2/)')
      write(iw,'(4x,a)') '========================'
      write(iw,'(4x,a)') 'RESP charges calculation'
      write(iw,'(4x,a)') '========================'
    case default
      error stop 'Unknown type of ESP charges calculation'
    end select

    call flush(iw)

!   Load basis set
    basis => infos%basis
    basis%atoms => infos%atoms

!   Allocate memory
    nbf = basis%nbf
    nbf2 = nbf*(nbf+1)/2
    nat = ubound(infos%atoms%zn, 1)
    npt = nat * sum(npt_layer)
    urohf = infos%control%scftype == 2 .or. infos%control%scftype == 3

    allocate(xyz(npt,3), &
             wt(npt), &
             pot(npt), &
             leb(maxval(npt_layer),3), &
             lebw(maxval(npt_layer)), &
             vdwrad(nat), &
             chg(nat), &
             neigh(nat), &
             den(nbf2), &
             stat=ok)
    if (ok /= 0) call show_message('Cannot allocate memory', WITH_ABORT)

!   Set up the grid
    nptcur = 0  ! Current number of point in a grid

!   Loop over layers and add spherical grid points on each atom to the molecular grid
    do layer = 1, nlayers

!     Get the atomic radii on which to place new grid layer
      vdwrad = ELEMENTS_VDW_RADII(int(infos%atoms%zn))*layers(layer)

!     Get grid
      nleb = npt_layer(layer)
      leb = 0
      lebw = 0
      call lebedev_get_grid(nleb, leb, lebw, typ_layer(layer))

!     Add new grid layer for each atom, remove inner points
      do i = 1, nat
        call add_atom_grid( &
          x=xyz(nptcur+1:,1), &
          y=xyz(nptcur+1:,2), &
          z=xyz(nptcur+1:,3), &
          wts=wt(nptcur+1:), &
          nadd=nadd, &
          atpts=leb(:nleb,:), &
          atwts=lebw(:nleb), &
          atoms_xyz=infos%atoms%xyz, &
          atoms_rad=vdwrad, &
          cur_atom=i, &
          neighbours=neigh)
        nptcur = nptcur + nadd
      end do
    end do

!   Set all grid weights to be 1 for now
    wt = 1

    deallocate(leb, lebw, vdwrad, neigh)

!   Get density
    if (urohf) then
      call data_has_tags(infos%dat, tags_beta, module_name, subroutine_name, WITH_ABORT)
      call tagarray_get_data(infos%dat, OQP_DM_A, dmat_a)
      call tagarray_get_data(infos%dat, OQP_DM_B, dmat_b)
      den = dmat_a + dmat_b
    else
      call data_has_tags(infos%dat, tags_alpha, module_name, subroutine_name, WITH_ABORT)
      call tagarray_get_data(infos%dat, OQP_DM_A, dmat_a)
      den = dmat_a
    end if

!   Compute electrostatic potential
    pot = 0
    call electrostatic_potential(basis, &
            x=xyz(:nptcur,1), &
            y=xyz(:nptcur,2), &
            z=xyz(:nptcur,3), &
            wt=wt(:nptcur), &
            d=den, &
            pot=pot)

    deallocate(den)

!   Add nuclei contribution to the potential
    call nuc_pot( &
            x=xyz(:nptcur,1), &
            y=xyz(:nptcur,2), &
            z=xyz(:nptcur,3), &
            w=wt(:nptcur), &
            at=infos%atoms%xyz, &
            q=infos%atoms%zn, &
            pot=pot(:nptcur))

!   Fit ESP charges
    call chg_fit_mk( &
            x=xyz(:nptcur,1), &
            y=xyz(:nptcur,2), &
            z=xyz(:nptcur,3), &
            w=wt(:nptcur), &
            at=infos%atoms%xyz, &
            pot=pot(:nptcur), &
            chgtot=real(infos%mol_prop%charge, dp), &
            chg=chg, &
            resp=restr, &
            q0=q0, &
            alpha=alpha &
    )

!   Print ESP charges
    call print_charges(infos, chg)

!   Compute RMS error of the ponential induced by ESP charges
    call check_charges( &
            x=xyz(:nptcur,1), &
            y=xyz(:nptcur,2), &
            z=xyz(:nptcur,3), &
            w=wt(:nptcur), &
            at=infos%atoms%xyz, &
            chg=chg, &
            pot=pot(:nptcur), &
            rms=rms)
    write(*,'(x,a,es20.6)') 'rms.err.=', rms

    close(iw)

  end subroutine oqp_resp_charges

!--------------------------------------------------------------------------------

!> @brief Add atomic-centered spherical grid to the molecular grid for ESP calculations
!> @param[in,out]  x            X coordinates of a molecular grid
!> @param[in,out]  y            Y coordinates of a molecular grid
!> @param[in,out]  z            Z coordinates of a molecular grid
!> @param[in,out]  wts          molecular grid weights
!> @param[out]     nadd         number of points added
!> @param[in]      atpts        atomic grid points (npts,3)
!> @param[in]      atwts        atomic grid weights
!> @param[in]      atoms_xyz    atomic coordinates
!> @param[in]      atoms_rad    atomic radii == radii of the spherical grid
!> @param[in]      cur_atom     id of the current atom
!> @param[in,out]  neighbours   temporary array to store current atom neighbours list
!
!> @detail Atoms are supposed to be speres with radii specified in atoms_rad.
!>         Neighbours are atoms, intesecting with the current atom, i.e. D_ij < R_i + R_j
!>         This subroutine adds only those points of the current atom, which are outside of spehres
!>         of all other atoms
!
!> @author   Vladimir Mironov
!
!     REVISION HISTORY:
!> @date _Mar, 2023_ Initial release
  subroutine add_atom_grid(x, y, z, wts, nadd, atpts, atwts, atoms_xyz, atoms_rad, cur_atom, neighbours)
    use precision, only: dp
    real(kind=dp), intent(inout) :: x(:), y(:), z(:), wts(:)
    real(kind=dp), intent(in) :: atpts(:,:), atwts(:)
    real(kind=dp), intent(in) :: atoms_xyz(:,:), atoms_rad(:)
    integer, intent(in) :: cur_atom
    integer, intent(inout) :: neighbours(:)
    integer, intent(out) :: nadd

    integer :: nat, natpts, nngh, i, j, nb
    real(kind=dp) :: cur_xyz(3), ptxyz(3), rcur, dist2
    logical :: add

    nadd = 0

    nat = ubound(atoms_xyz, 2)
    natpts = ubound(atpts, 1)
    cur_xyz = atoms_xyz(:,cur_atom)
    rcur = atoms_rad(cur_atom)

    ! Find all neighbours of the current atom: D_ij < R^{vdw}_i + R^{vdw}_j
    nngh = 0
    do i = 1, nat
      if (i==cur_atom) cycle
      dist2 = sum((atoms_xyz(:,i) - cur_xyz(:))**2) - (atoms_rad(i)+rcur)**2
      if (dist2 < 0) then
        nngh = nngh + 1
        neighbours(nngh) = i
      end if
    end do

    ! Add only the points which are outside of the VDW shell of all other atoms
    do j = 1, natpts
      add = .true.
      ptxyz = cur_xyz + rcur*atpts(j,:)
      ! Check if grid is outside all other atoms VDW spheres
      do i = 1, nngh
        nb = neighbours(i)
        dist2 = sum((atoms_xyz(:,nb) - ptxyz)**2)
        add = dist2 > atoms_rad(nb)**2
        if (.not.add) exit
      end do

      ! Add new point to the grid, increase the counter
      if (add) then
        nadd = nadd + 1
        x(nadd) = ptxyz(1)
        y(nadd) = ptxyz(2)
        z(nadd) = ptxyz(3)
        wts(nadd) = atwts(j)
      end if
    end do

  end subroutine add_atom_grid

!--------------------------------------------------------------------------------

!> @brief Compute nuclear potential of atoms on a grid
!> @param[in]      x            X coordinates of a molecular grid
!> @param[in]      y            Y coordinates of a molecular grid
!> @param[in]      z            Z coordinates of a molecular grid
!> @param[in]      w            molecular grid weights
!> @param[in]      at           coordinates of atoms
!> @param[in]      q            atomic nuclei charges
!> @param[out]     pot          values of nuclear potential on a grid
!
!> @author   Vladimir Mironov
!
!     REVISION HISTORY:
!> @date _Mar, 2023_ Initial release
  subroutine nuc_pot(x, y, z, w, at, q, pot)
    use precision, only: dp
    real(kind=dp), intent(in) :: x(:), y(:), z(:), w(:), at(:,:), q(:)
    real(kind=dp), intent(out) :: pot(:)

    integer :: i, j, nat, npts
    nat = ubound(q, 1)
    npts = ubound(x, 1)

    ! i - atoms
    ! j - grid pts
    ! pot_j = \sum_i^{N_{at}} w_j * Q_i / D_ji
    do j = 1, npts
      do i = 1, nat
        pot(j) = pot(j) &
               + w(j)*q(i)/norm2(at(:,i) - [x(j), y(j), z(j)])
      end do
    end do

  end subroutine nuc_pot

!--------------------------------------------------------------------------------

!> @brief Compute ESP charges using generalized least-squares fit in LAPACK
!> @param[in]      x            X coordinates of a molecular grid
!> @param[in]      y            Y coordinates of a molecular grid
!> @param[in]      z            Z coordinates of a molecular grid
!> @param[in]      w            molecular grid weights (not used for now)
!> @param[in]      at           coordinates of atoms
!> @param[in]      pot          values of nuclear potential on a grid
!> @param[in]      chgtot       total charge of the system
!> @param[out]     chg          ESP charges of atoms
!
!> @author   Vladimir Mironov
!
!     REVISION HISTORY:
!> @date _Mar, 2023_ Initial release
  subroutine chg_fit_glsq(x, y, z, w, at, pot, chgtot, chg)
    use precision, only: dp
    real(kind=dp), intent(in) :: x(:), y(:), z(:), w(:), at(:,:), pot(:), chgtot
    real(kind=dp), intent(inout) :: chg(:)

    integer :: i, j, nat, npts
    real(kind=dp), allocatable :: a(:,:), b(:), c(:)
    real(kind=dp), allocatable :: work(:)
    integer :: lwork, info
    real(kind=dp) :: rwork(1), d(1)

    nat = ubound(at, 2)
    npts = ubound(x, 1)

    allocate(b(nat), c(npts), a(npts,nat), source=0.0d0)

    ! Solve the problem:
    !    minimize || pot - D*chg ||_2
    !    subject to \sum(chg) = chgtot
    ! In LAPACK we need to formulate the contraint in the form B*x = d
    ! which will be:
    !   (1,1,1...1)^T * chg = chgtot

    ! Compute the matrix of the inverse distances
    do i = 1, nat
      do j = 1, npts
        a(j,i) = 1/norm2(at(:,i) - [x(j), y(j), z(j)])
      end do
    end do

    ! Copy target values for LAPACK
    c(:npts) = pot

    ! Set up constraints
    b = 1
    d(1) = chgtot

    ! Query required workspace
    call dgglse(npts, nat, 1, a, npts, b, 1, c, d, chg, rwork, -1, info)

    lwork = nint(rwork(1),8)
    allocate(work(lwork))

    ! Solve the problem
    call dgglse(npts, nat, 1, a, npts, b, 1, c, d, chg, work, lwork, info)

  end subroutine chg_fit_glsq

!--------------------------------------------------------------------------------

!> @brief Compute ESP charges using Merz-Kollman algorithm
!> @details This subroutine uses the Lagrangian proposed by MK
!>   See the following paper for details:
!>   B.Besler, K. Merz, P.Kollman, J. Comput. Chem., 11, 4, 431-439 (1990)
!> @param[in]      x            X coordinates of a molecular grid
!> @param[in]      y            Y coordinates of a molecular grid
!> @param[in]      z            Z coordinates of a molecular grid
!> @param[in]      w            molecular grid weights (not used for now)
!> @param[in]      at           coordinates of atoms
!> @param[in]      pot          values of nuclear potential on a grid
!> @param[in]      chgtot       total charge of the system
!> @param[out]     chg          ESP charges of atoms
!> @param[in]      resp         if .true., run restrainted ESP (RESP) fit
!> @param[in]      q0           RESP target charges
!> @param[in]      alpha        RESP constraints
!
!> @author   Vladimir Mironov
!
!     REVISION HISTORY:
!> @date _Mar, 2023_ Initial release
  subroutine chg_fit_mk(x, y, z, w, at, pot, chgtot, chg, resp, q0, alpha)
    use precision, only: dp
    real(kind=dp), intent(in) :: x(:), y(:), z(:), w(:), at(:,:), pot(:), chgtot
    real(kind=dp), intent(inout) :: chg(:)
    logical, optional, intent(in) :: resp
    real(kind=dp), optional, intent(in) :: q0(:), alpha(:)

    integer, parameter :: RESP_MAX_ITER = 30
    real(kind=dp), parameter :: RESP_TOL = 1.0d-5
    logical, parameter :: debug = .false.

    integer :: i, j, nat, npts
    real(kind=dp), allocatable :: a(:,:), b(:), r(:,:)
    real(kind=dp), allocatable :: a_bak(:,:), b_bak(:)
    integer, allocatable :: ipiv(:)
    integer :: info
    logical :: do_resp
    real(kind=dp) :: diff

    do_resp = .false.
    if (present(resp)) do_resp = resp

    nat = ubound(at, 2)
    npts = ubound(x, 1)

    allocate(a(nat+1,nat+1), b(nat+1), r(nat,npts), ipiv(nat+1))

    ! Compute the matrix of the inverse distances
    do j = 1, npts
      do i = 1, nat
        r(i,j) = 1/norm2(at(:,i) - [x(j), y(j), z(j)])
      end do
    end do

    ! Compute the main part of the Lagrangian
    call dgemm('n','t', &
            nat, nat, npts, &
            1.0_dp, r, nat, &
                    r, nat, &
            0.0_dp, a, nat+1)

    ! Compute the part of the Lagrangian corresponding to the constraint:
    !   sum(chg) = chgtot
    a(:,nat+1) = 1
    a(nat+1,:) = 1
    a(nat+1,nat+1) = 0
    ! Compute the vector b
    call dgemv('n', nat, npts, 1.0_dp, r, nat, pot, 1, 0.0_dp, b, 1)

    b(nat+1) = chgtot

    ! We don't need the inverse distance matrix anymore
    deallocate(r)

    ! Save the original A and b in case of RESP, because
    ! GESV will destroy them
    if (do_resp) then
      a_bak = a
      b_bak = b
    end if

    ! Solve linear equation A*chg = b
    call dgesv(nat+1, 1, a, nat+1, ipiv, b, nat+1, info)

    chg(:nat) = b(:nat)

    ! Iterative solution of RESP equations
    if (do_resp) then
      do j = 1, RESP_MAX_ITER
        ! Recover original A and b
        a = a_bak
        b = b_bak

        ! Add harmonic restraint contribution:
        do i = 1, nat
          a(i,i) = a(i,i) + 2*alpha(i)*(q0(i)-chg(i))
        end do
        b(:nat) = b(:nat) + 2*q0(:nat)*alpha(:nat)*(q0(:nat)-chg(:nat))

        ! Solve the equation:
        call dgesv(nat+1, 1, a, nat+1, ipiv, b, nat+1, info)

        ! Check convergence of charges
        diff = maxval(abs(chg(:nat)-b(:nat)))
        if (debug) write(*,'(X,A10,I4,A10,ES10.3)') 'resp iter=', j, ' err=', diff

        ! Copy the solution, exit if converged
        chg(:nat) = b(:nat)
        if (diff < RESP_TOL) exit

      end do
    end if

  end subroutine chg_fit_mk

!--------------------------------------------------------------------------------

!> @brief Compute an RMS error of the pontential induced by atomic charges
!>
!> @param[in]      x            X coordinates of a molecular grid
!> @param[in]      y            Y coordinates of a molecular grid
!> @param[in]      z            Z coordinates of a molecular grid
!> @param[in]      w            molecular grid weights
!> @param[in]      at           coordinates of atoms
!> @param[in]      chg          atomic partial charges
!> @param[in]      pot          reference ponential values on a grid
!> @param[out]     rms          RMS error
!
!> @author   Vladimir Mironov
!
!     REVISION HISTORY:
!> @date _Mar, 2023_ Initial release
  subroutine check_charges(x, y, z, w, at, chg, pot, rms)
    use precision, only: dp
    real(kind=dp), intent(in)  :: x(:), y(:), z(:), w(:), at(:,:), chg(:)
    real(kind=dp), intent(in)  :: pot(:)
    real(kind=dp), intent(out) :: rms

    integer :: i, j, nat, npts
    real(kind=dp) :: v, vtot

    npts = ubound(x, 1)
    nat = ubound(chg, 1)

    vtot = 0
    do j = 1, npts
      v = 0
      do i = 1, nat
        v = v &
          + w(j)*chg(i)/norm2(at(:,i) - [x(j), y(j), z(j)])
      end do
      vtot = vtot + (v-pot(j))**2
    end do

    rms = sqrt(vtot/npts)

  end subroutine check_charges

!--------------------------------------------------------------------------------

!> @brief Print partial charges
!>
!> @param[in]      infos        OQP handle
!> @param[in]      chg          atomic partial charges
!
!> @author   Vladimir Mironov
!
!     REVISION HISTORY:
!> @date _Mar, 2023_ Initial release
  subroutine print_charges(infos, chg)
    use precision, only: dp
    use elements, only: ELEMENTS_SHORT_NAME
    use types, only: information
    type(information), intent(in) :: infos
    real(kind=dp), intent(in) :: chg(:)

    integer :: i, elem, nat
    nat = ubound(infos%atoms%zn, 1)

    write(*,'(/,30("^"))')
    write(*,'(/a8,a8,a14)') '#', 'Name', 'Charge'
    write(*,'(30("-"))')

    do i = 1, nat
      elem = nint(infos%atoms%zn(i))
      write(*,'(i8,a8,f14.6)') i, ELEMENTS_SHORT_NAME(elem), chg(i)
    end do

    write(*,'(30("="))')

  end subroutine print_charges

end module resp_mod
