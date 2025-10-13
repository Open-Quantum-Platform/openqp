module tdhf_mrsf_z_vector_mod

  implicit none

  character(len=*), parameter :: module_name = "tdhf_mrsf_z_vector_mod"

  ! Module-level work arrays for GMRES to avoid repeated allocation
  real(kind=8), allocatable :: gmres_wrk1(:,:), gmres_wrk2(:,:), gmres_wrk3(:,:)
  real(kind=8), allocatable, target :: gmres_pa(:,:,:)
  real(kind=8), allocatable :: gmres_ab1_mo_a(:,:), gmres_ab1_mo_b(:,:)
  logical :: gmres_work_allocated = .false.
  integer :: gmres_nbf = 0
  integer :: gmres_nocca = 0
  integer :: gmres_noccb = 0

contains

  ! Initialize GMRES work arrays
  subroutine init_gmres_work(nbf, nocca, noccb)
    use messages, only: show_message, with_abort
    implicit none
    integer, intent(in) :: nbf, nocca, noccb
    integer :: nvira, nvirb, ok
    
    if (gmres_work_allocated) then
      ! Check if dimensions match
      if (gmres_nbf == nbf .and. gmres_nocca == nocca .and. gmres_noccb == noccb) then
        return  ! Arrays already allocated with correct size
      else
        ! Deallocate old arrays before reallocating
        call cleanup_gmres_work()
      end if
    end if
    
    nvira = nbf - nocca
    nvirb = nbf - noccb
    
    allocate(gmres_wrk1(nbf,nbf), &
             gmres_wrk2(nbf,nbf), &
             gmres_wrk3(nbf,nbf), &
             gmres_pa(nbf,nbf,2), &
             gmres_ab1_mo_a(nocca,nvira), &
             gmres_ab1_mo_b(noccb,nvirb), &
             stat=ok)
    
    if (ok /= 0) then
      call show_message('Cannot allocate GMRES work arrays', with_abort)
    end if
    
    gmres_work_allocated = .true.
    gmres_nbf = nbf
    gmres_nocca = nocca
    gmres_noccb = noccb
    
  end subroutine init_gmres_work
  
  ! Cleanup GMRES work arrays
  subroutine cleanup_gmres_work()
    implicit none
    
    if (allocated(gmres_wrk1)) deallocate(gmres_wrk1)
    if (allocated(gmres_wrk2)) deallocate(gmres_wrk2)
    if (allocated(gmres_wrk3)) deallocate(gmres_wrk3)
    if (allocated(gmres_pa)) deallocate(gmres_pa)
    if (allocated(gmres_ab1_mo_a)) deallocate(gmres_ab1_mo_a)
    if (allocated(gmres_ab1_mo_b)) deallocate(gmres_ab1_mo_b)
    
    gmres_work_allocated = .false.
    gmres_nbf = 0
    gmres_nocca = 0
    gmres_noccb = 0
    
  end subroutine cleanup_gmres_work

  ! GMRES solver for the z-vector equation
  subroutine gmres_solve(apply_operator, apply_precond, b, x, n, restart, max_iter, tol, &
                         infos, basis, molGrid, int2_driver, &
                         nocca, noccb, nbf, mo_a, mo_b, mo_energy_a, &
                         fa, fb, scale_exch, dft, error_out, iter_out, iw)
    use precision, only: dp
    use types, only: information
    use basis_tools, only: basis_set
    use int2_compute, only: int2_compute_t
    use mod_dft_molgrid, only: dft_grid_t
    
    implicit none
    
    interface
      subroutine apply_operator(x_in, x_out, infos, basis, molGrid, int2_driver, &
                               nocca, noccb, nbf, mo_a, mo_b, mo_energy_a, &
                               fa, fb, scale_exch, dft)
        use precision, only: dp
        use types, only: information
        use basis_tools, only: basis_set
        use int2_compute, only: int2_compute_t
        use mod_dft_molgrid, only: dft_grid_t
        real(kind=dp), intent(in) :: x_in(:)
        real(kind=dp), intent(out) :: x_out(:)
        type(information), intent(inout) :: infos
        type(basis_set), pointer :: basis
        type(dft_grid_t), intent(inout) :: molGrid
        type(int2_compute_t), intent(inout) :: int2_driver
        integer, intent(in) :: nocca, noccb, nbf
        real(kind=dp), intent(in) :: mo_a(:,:), mo_b(:,:), mo_energy_a(:)
        real(kind=dp), intent(in) :: fa(:,:), fb(:,:), scale_exch
        logical, intent(in) :: dft
      end subroutine
      
      subroutine apply_precond(x_in, x_out)
        use precision, only: dp
        real(kind=dp), intent(in) :: x_in(:)
        real(kind=dp), intent(out) :: x_out(:)
      end subroutine
    end interface
    
    real(kind=dp), intent(in) :: b(:)
    real(kind=dp), intent(inout) :: x(:)
    integer, intent(in) :: n, restart, max_iter, iw
    real(kind=dp), intent(in) :: tol
    type(information), intent(inout) :: infos
    type(basis_set), pointer :: basis
    type(dft_grid_t), intent(inout) :: molGrid
    type(int2_compute_t), intent(inout) :: int2_driver
    integer, intent(in) :: nocca, noccb, nbf
    real(kind=dp), intent(in) :: mo_a(:,:), mo_b(:,:), mo_energy_a(:)
    real(kind=dp), intent(in) :: fa(:,:), fb(:,:), scale_exch
    logical, intent(in) :: dft
    real(kind=dp), intent(out) :: error_out
    integer, intent(out) :: iter_out
    
    ! Local variables
    real(kind=dp), allocatable :: V(:,:)     ! Krylov basis
    real(kind=dp), allocatable :: H(:,:)     ! Hessenberg matrix
    real(kind=dp), allocatable :: c(:), s(:) ! Givens rotation coefficients
    real(kind=dp), allocatable :: g(:)       ! RHS for least squares
    real(kind=dp), allocatable :: y(:)       ! Solution of least squares
    real(kind=dp), allocatable :: r(:)       ! Residual
    real(kind=dp), allocatable :: w(:)       ! Work vector
    real(kind=dp), allocatable :: Ax(:)      ! A*x
    
    real(kind=dp) :: beta, h_ij, temp, error, error_initial
    integer :: i, j, k, iter, m, restart_count, inner_iter
    logical :: converged
    
    ! Initialize GMRES work arrays ONCE at the beginning
    call init_gmres_work(nbf, nocca, noccb)
    
    ! Allocate workspace
    m = min(restart, n)
    allocate(V(n, m+1))
    allocate(H(m+1, m))
    allocate(c(m))
    allocate(s(m))
    allocate(g(m+1))
    allocate(y(m))
    allocate(r(n))
    allocate(w(n))
    allocate(Ax(n))
    
    iter_out = 0
    restart_count = 0
    converged = .false.
    
    write(iw,'(/," GMRES Solver Parameters:")')
    write(iw,'("   Problem size        : ", I8)') n
    write(iw,'("   Restart dimension   : ", I4)') m
    write(iw,'("   Max iterations      : ", I4)') max_iter
    write(iw,'("   Convergence tol     : ", 1p,e10.3)') tol
    write(iw,'(/," Iteration   Inner  Residual Norm   Reduction")')
    write(iw,'(" ---------   -----  -------------   ---------")')
    call flush(iw)
    
    ! Compute initial residual
    call apply_operator(x, Ax, infos, basis, molGrid, int2_driver, &
                       nocca, noccb, nbf, mo_a, mo_b, mo_energy_a, &
                       fa, fb, scale_exch, dft)
    r = b - Ax
    error_initial = sqrt(dot_product(r, r))
    
    ! Outer iteration loop (restarts)
    do iter = 1, max_iter
      
      restart_count = restart_count + 1
      
      ! Compute initial residual r = b - A*x
      call apply_operator(x, Ax, infos, basis, molGrid, int2_driver, &
                         nocca, noccb, nbf, mo_a, mo_b, mo_energy_a, &
                         fa, fb, scale_exch, dft)
      r = b - Ax
      
      ! Apply preconditioner to residual
      call apply_precond(r, V(:,1))
      
      beta = sqrt(dot_product(V(:,1), V(:,1)))
      
      ! Check for convergence
      error = beta
      if (iter == 1) then
        write(iw,'(I6,8x,"  0",2x,1p,F13.8,1x,F13.8)') &
              restart_count, error, error/error_initial
      end if
      
      if (error < tol) then
        error_out = error
        iter_out = iter
        converged = .true.
        exit
      end if
      
      V(:,1) = V(:,1) / beta
      g = 0.0_dp
      g(1) = beta
      
      ! Reset H matrix for this restart
      H = 0.0_dp
      
      ! Arnoldi process
      inner_iter = 0
      do j = 1, m
        inner_iter = j
        
        ! Apply operator to V_j
        call apply_operator(V(:,j), w, infos, basis, molGrid, int2_driver, &
                           nocca, noccb, nbf, mo_a, mo_b, mo_energy_a, &
                           fa, fb, scale_exch, dft)
        
        ! Apply preconditioner
        call apply_precond(w, V(:,j+1))
        
        ! Modified Gram-Schmidt orthogonalization
        do i = 1, j
          H(i,j) = dot_product(V(:,j+1), V(:,i))
          V(:,j+1) = V(:,j+1) - H(i,j) * V(:,i)
        end do
        
        ! Reorthogonalization for numerical stability
        do i = 1, j
          temp = dot_product(V(:,j+1), V(:,i))
          H(i,j) = H(i,j) + temp
          V(:,j+1) = V(:,j+1) - temp * V(:,i)
        end do
        
        H(j+1,j) = sqrt(dot_product(V(:,j+1), V(:,j+1)))
        
        ! Check for breakdown
        if (abs(H(j+1,j)) < 1.0d-12) then
          write(iw,'(" GMRES: Lucky breakdown at iteration ", I3)') j
          inner_iter = j
          exit
        end if
        
        V(:,j+1) = V(:,j+1) / H(j+1,j)
        
        ! Apply previous Givens rotations
        do i = 1, j-1
          temp = c(i) * H(i,j) + s(i) * H(i+1,j)
          H(i+1,j) = -s(i) * H(i,j) + c(i) * H(i+1,j)
          H(i,j) = temp
        end do
        
        ! Compute new Givens rotation
        call givens_rotation(H(j,j), H(j+1,j), c(j), s(j))
        
        ! Apply new Givens rotation
        H(j,j) = c(j) * H(j,j) + s(j) * H(j+1,j)
        H(j+1,j) = 0.0_dp
        
        temp = c(j) * g(j) + s(j) * g(j+1)
        g(j+1) = -s(j) * g(j) + c(j) * g(j+1)
        g(j) = temp
        
        ! Check convergence
        error = abs(g(j+1)) 
        iter_out = iter_out + 1
        
        ! Print progress every 5 inner iterations or at convergence
        if (mod(j, 5) == 0 .or. error < tol .or. j == m) then
          write(iw,'(I6,8x,I3,2x,1p,F13.8,1x,F13.8)') &
                restart_count, j, error, error/error_initial
          call flush(iw)
        end if
        
        if (error < tol) then
          converged = .true.
          inner_iter = j
          exit
        end if
        
        if (iter_out >= max_iter) then
          inner_iter = j
          exit
        end if
      end do
      
      ! Solve upper triangular system for y
      call back_substitution(H(1:inner_iter,1:inner_iter), g(1:inner_iter), &
                             y(1:inner_iter), inner_iter)
      
      ! Update solution: x = x + V*y
      do i = 1, inner_iter
        x = x + y(i) * V(:,i)
      end do
      
      error_out = error
      
      if (converged .or. iter_out >= max_iter) exit
      
      ! Print restart information
      if (.not. converged .and. inner_iter == m) then
        write(iw,'(" GMRES: Restarting (restart #", I3, ")")') restart_count
        call flush(iw)
      end if
      
    end do
    
    ! Final status
    write(iw,'(" ---------   -----  -------------   ---------")')
    if (converged) then
      write(iw,'(" GMRES converged in ", I4, " iterations (", I3, " restarts)")') &
            iter_out, restart_count-1
      write(iw,'(" Final residual norm: ", 1p,e13.6)') error_out
    else
      write(iw,'(" GMRES did not converge within ", I4, " iterations")') max_iter
      write(iw,'(" Final residual norm: ", 1p,e13.6)') error_out
    end if
    write(iw,'(" Relative reduction : ", 1p,e13.6)') error_out/error_initial
    call flush(iw)
    
    ! Clean up local arrays
    deallocate(V, H, c, s, g, y, r, w, Ax)
    
    ! NOTE: Do NOT clean up GMRES work arrays here - they will be cleaned in main routine
    
  contains
    
    subroutine givens_rotation(a, b, c, s)
      real(kind=dp), intent(in) :: a, b
      real(kind=dp), intent(out) :: c, s
      real(kind=dp) :: r
      
      if (abs(b) < 1.0d-14) then
        c = 1.0_dp
        s = 0.0_dp
      else
        r = sqrt(a*a + b*b)
        c = a / r
        s = b / r
      end if
    end subroutine givens_rotation
    
    subroutine back_substitution(A, b, x, n)
      integer, intent(in) :: n
      real(kind=dp), intent(in) :: A(n,n), b(n)
      real(kind=dp), intent(out) :: x(n)
      integer :: i, j
      
      x(n) = b(n) / A(n,n)
      do i = n-1, 1, -1
        x(i) = b(i)
        do j = i+1, n
          x(i) = x(i) - A(i,j) * x(j)
        end do
        x(i) = x(i) / A(i,i)
      end do
    end subroutine back_substitution
    
  end subroutine gmres_solve

  ! Apply the z-vector operator (A*x) - Modified to use module-level arrays
  subroutine apply_z_operator(x_in, x_out, infos, basis, molGrid, int2_driver, &
                              nocca, noccb, nbf, mo_a, mo_b, mo_energy_a, &
                              fa, fb, scale_exch, dft)
    use precision, only: dp
    use types, only: information
    use basis_tools, only: basis_set
    use int2_compute, only: int2_compute_t
    use tdhf_lib, only: int2_tdgrd_data_t
    use tdhf_sf_lib, only: sfrogen, sfrolhs
    use mod_dft_gridint_fxc, only: utddft_fxc
    use mathlib, only: symmetrize_matrix, orthogonal_transform
    use mod_dft_molgrid, only: dft_grid_t
    use tdhf_lib, only: mntoia
    
    implicit none
    
    real(kind=dp), intent(in) :: x_in(:)
    real(kind=dp), intent(out) :: x_out(:)
    type(information), intent(inout) :: infos
    type(basis_set), pointer :: basis
    type(dft_grid_t), intent(inout) :: molGrid
    type(int2_compute_t), intent(inout) :: int2_driver
    integer, intent(in) :: nocca, noccb, nbf
    real(kind=dp), intent(in) :: mo_a(:,:), mo_b(:,:), mo_energy_a(:)
    real(kind=dp), intent(in) :: fa(:,:), fb(:,:), scale_exch
    logical, intent(in) :: dft
    
    ! Local variables
    real(kind=dp), pointer :: ab1(:,:,:)
    type(int2_tdgrd_data_t), allocatable, target :: int2_data
    integer :: nvira, nvirb
    
    nvira = nbf - nocca
    nvirb = nbf - noccb
    
    ! Ensure work arrays are initialized and have correct dimensions
    if (.not. gmres_work_allocated .or. &
        gmres_nbf /= nbf .or. &
        gmres_nocca /= nocca .or. &
        gmres_noccb /= noccb) then
      call init_gmres_work(nbf, nocca, noccb)
    end if
    
    ! Clear work arrays
    gmres_wrk1 = 0.0_dp
    gmres_wrk2 = 0.0_dp
    gmres_wrk3 = 0.0_dp
    gmres_pa = 0.0_dp
    gmres_ab1_mo_a = 0.0_dp
    gmres_ab1_mo_b = 0.0_dp
    
    ! Generate density matrices from x_in
    call sfrogen(gmres_wrk1, gmres_wrk2, x_in, nocca, noccb)
    
    ! Transform to AO basis
    call orthogonal_transform('t', nbf, mo_a, gmres_wrk1, gmres_pa(:,:,1), gmres_wrk3)
    call orthogonal_transform('t', nbf, mo_b, gmres_wrk2, gmres_pa(:,:,2), gmres_wrk3)
    
    ! Initialize ERI calculation with proper allocation
    allocate(int2_data)
    int2_data = int2_tdgrd_data_t( &
        d2 = gmres_pa, &
        int_apb = .true., &
        int_amb = .false., &
        tamm_dancoff = .false., &
        scale_exchange = scale_exch)
    
    call int2_driver%run(int2_data, &
          cam=dft.and.infos%dft%cam_flag, &
          alpha=infos%dft%cam_alpha, &
          beta=infos%dft%cam_beta,&
          mu=infos%dft%cam_mu)
    ab1 => int2_data%apb(:,:,:,1)
    
    call symmetrize_matrix(gmres_pa(:,:,1), nbf)
    call symmetrize_matrix(gmres_pa(:,:,2), nbf)
    
    if (dft) then
      call utddft_fxc( &
          basis = basis, &
          molGrid = molGrid, &
          isVecs = .true., &
          wfa = mo_a, &
          wfb = mo_b, &
          fxa = ab1(:,:,1:1), &
          fxb = ab1(:,:,2:2), &
          dxa = gmres_pa(:,:,1:1), &
          dxb = gmres_pa(:,:,2:2), &
          nmtx = 1, &
          threshold = 1.0d-15, &
          infos = infos)
    end if
    
    ! Transform to MO basis - Fixed to use correct mo_b for beta
    call mntoia(ab1(:,:,1), gmres_ab1_mo_a, mo_a, mo_a, nocca, nocca)
    call mntoia(ab1(:,:,2), gmres_ab1_mo_b, mo_b, mo_b, noccb, noccb)
    
    ! Apply the operator
    call sfrolhs(x_out, x_in, mo_energy_a, fa, fb, gmres_ab1_mo_a, gmres_ab1_mo_b, &
                 nocca, noccb)
    
    call int2_data%clean()
    deallocate(int2_data)
    
  end subroutine apply_z_operator

  ! Apply preconditioner (simple diagonal preconditioner)
  subroutine apply_z_precond(x_in, x_out, xminv)
    use precision, only: dp
    implicit none
    real(kind=dp), intent(in) :: x_in(:), xminv(:)
    real(kind=dp), intent(out) :: x_out(:)
    
    x_out = xminv * x_in
    
  end subroutine apply_z_precond

  subroutine tdhf_mrsf_z_vector_C(c_handle) bind(C, name="tdhf_mrsf_z_vector")
    use c_interop, only: oqp_handle_t, oqp_handle_get_info
    use types, only: information
    type(oqp_handle_t) :: c_handle
    type(information), pointer :: inf
    inf => oqp_handle_get_info(c_handle)
    call tdhf_mrsf_z_vector(inf)
  end subroutine tdhf_mrsf_z_vector_C

  subroutine tdhf_mrsf_z_vector(infos)
    use precision, only: dp
    use io_constants, only: iw
    use oqp_tagarray_driver

    use types, only: information
    use strings, only: Cstring, fstring
    use basis_tools, only: basis_set
    use messages, only: show_message, with_abort
    use util, only: measure_time

    use int2_compute, only: int2_compute_t
    use tdhf_lib, only: int2_td_data_t
    use tdhf_lib, only: int2_tdgrd_data_t
    use tdhf_mrsf_lib, only: int2_mrsf_data_t
    use tdhf_lib, only: iatogen, mntoia
    use tdhf_sf_lib, only: sfrorhs, &
      sfromcal, sfrogen, sfrolhs, pcgrbpini, &
      pcgb, sfropcal, sfdmat
    use dft, only: dft_initialize, dftclean
    use mod_dft_gridint_fxc, only: utddft_fxc
    use mathlib, only: symmetrize_matrix, orthogonal_transform, &
            orthogonal_transform_sym
    use mod_dft_molgrid, only: dft_grid_t
    use mathlib, only: pack_matrix, unpack_matrix

    use tdhf_mrsf_lib, only: &
      mrinivec, mrsfcbc, mrsfxvec, mrsfsp, mrsfrowcal, &
      mrsfqrorhs, mrsfqropcal, mrsfqrowcal
    use oqp_linalg
    use printing, only: print_module_info


    implicit none

    character(len=*), parameter :: subroutine_name = "tdhf_mrsf_z_vector"

    type(basis_set), pointer :: basis
    type(information), target, intent(inout) :: infos

    integer :: ok

    real(kind=dp), allocatable :: ab1_mo_a(:,:)
    real(kind=dp), allocatable :: ab1_mo_b(:,:)
    real(kind=dp), allocatable :: xm(:)
    real(kind=dp), pointer :: ab2(:,:,:)
    real(kind=dp), pointer :: ab1(:,:,:)
    real(kind=dp), allocatable :: fa(:,:), fb(:,:)
    real(kind=dp), pointer :: bvec(:,:,:)
    real(kind=dp), pointer :: wmo(:,:)
    real(kind=dp), allocatable :: bvec_mo_d(:,:)
    real(kind=dp), allocatable, target :: &
      fmrst1(:,:,:,:)
    real(kind=dp), pointer :: fmrst2(:,:,:,:)

    integer :: nocca, nvira, noccb, nvirb
    integer :: nbf, nbf_tri
    integer :: iter, gmres_iter
    real(kind=dp) :: cnvtol, scale_exch, scale_exch2
    logical :: roref = .false.
    integer :: mrst

    type(int2_compute_t) :: int2_driver
    type(int2_mrsf_data_t), allocatable, target :: int2_data_st
    type(int2_td_data_t), allocatable, target :: int2_data_q
    class(int2_td_data_t), allocatable, target :: int2_data
    type(dft_grid_t) :: molGrid

  ! scr data
    real(kind=dp), allocatable, target :: wrk1(:,:), wrk2(:,:), wrk3(:,:)
    real(kind=dp), pointer :: wrk1t(:)

  ! SF-TD Gradient data
    real(kind=dp), allocatable :: &
      rhs(:), lhs(:), xminv(:), xk(:), pk(:), errv(:), &
      hxa(:,:), hxb(:,:), tij(:,:), ppija(:,:), ppijb(:,:), tab(:,:)
    real(kind=dp), allocatable, target :: pa(:,:,:)
    integer :: nsocc, lzdim, xvec_dim

  ! General data
    real(kind=dp) :: alpha, error
    character(len=10) :: solver_name

    logical :: dft
    integer :: scf_type, mol_mult, target_state

    ! tagarray
    real(kind=dp), contiguous, pointer :: &
      fock_a(:), mo_a(:,:), mo_energy_a(:), &
      fock_b(:), mo_b(:,:), &
      td_p(:,:), td_t(:,:), ta(:), tb(:), td_abxc(:,:), &
      td_mrsf_den(:,:,:), bvec_mo(:,:), wao(:), mrsf_energies(:)
    character(len=*), parameter :: tags_alloc(4) = (/ character(len=80) :: &
      OQP_WAO, OQP_td_mrsf_density, OQP_td_p, OQP_td_abxc /)
    character(len=*), parameter :: tags_required(8) = (/ character(len=80) :: &
      OQP_FOCK_A, OQP_E_MO_A, OQP_VEC_MO_A, OQP_FOCK_B, OQP_VEC_MO_B, OQP_td_bvec_mo, OQP_td_t, &
      OQP_td_energies /)

    mol_mult = infos%mol_prop%mult
    if (mol_mult/=3) call show_message(&
            'MRSF-TDDFT are available for ROHF/UHF ref.&
            &with ONLY triplet multiplicity(mult=3)', with_abort)

    scf_type = infos%control%scftype
    if (scf_type==3) roref = .true.

    dft = infos%control%hamilton == 20

  ! Files open
  ! 3. LOG: Write: Main output file
    open (unit=iw, file=infos%log_filename, position="append")
  !
    call print_module_info('MRSF_TDHF_Z_Vector','Solving Z-Vector for MRSF-TDDFT')

  ! Readings

  ! Load basis set
    basis => infos%basis
    basis%atoms => infos%atoms

    nbf = basis%nbf
    nbf_tri = nbf*(nbf+1)/2

    if (dft) call dft_initialize(infos, basis, molGrid)

  ! Parameter it should be inputed later
    mrst = infos%tddft%mult
    cnvtol = infos%tddft%zvconv

    nocca = infos%mol_prop%nelec_A
    nvira = nbf-noccA
    noccb = infos%mol_prop%nelec_B
    nvirb = nbf-noccb
    nsocc = nocca-noccb
    lzdim = noccb*(nsocc+nvira)+nsocc*nvira

    if(mrst==1 .or. mrst==3) then
      xvec_dim = nocca*nvirb
      allocate(&
    ! for Z-vector
        fmrst1(1,7,nbf,nbf), &
        bvec_mo_d(xvec_dim,1), &
        hxa(nbf,nocca), &
        hxb(nbf,nbf), &
    ! for gradient
        tij(nocca,nocca), &
        tab(nvirb,nvirb), &
        stat=ok, &
        source=0.0_dp)
    else if(mrst==5) then
      xvec_dim = noccb*nvira
      allocate(&
    ! for Z-vector
        hxa(nbf,nbf), &
        hxb(nbf,noccb), &
  ! for gradient
        tij(noccb,noccb), &
        tab(nvira,nvira), &
        stat=ok, &
        source=0.0_dp)
    endif
    if( ok/=0 ) call show_message('Cannot allocate memory', with_abort)

    allocate(&
  ! for Z-vector
      xminv(lzdim), &
      rhs(lzdim), &
      lhs(lzdim), &
      xm(lzdim), &
      xk(lzdim), &
      pk(lzdim), &
      errv(lzdim), &
   ! For gradient
      pa(nbf,nbf,2), &
      ppija(nocca,nocca), &
      ppijb(noccb,noccb), &
   ! Allocate TDDFT variables
      fa(nbf,nbf), &           ! Temporary matrix for diagonalization
      fb(nbf,nbf), &           ! Temporary matrix for diagonalization
      ab1_mo_a(nocca,nvira), &
      ab1_mo_b(noccb,nvirb), &
!   For scratch
      wrk1(nbf,nbf), &
      wrk2(nbf,nbf), &
      wrk3(nbf,nbf), &
      stat=ok, &
      source=0.0_dp)

    if( ok/=0 ) call show_message('Cannot allocate memory', with_abort)

    call infos%dat%remove_records(tags_alloc)

    call infos%dat%reserve_data(OQP_WAO, TA_TYPE_REAL64, nbf_tri, comment=OQP_WAO_comment)
    call infos%dat%reserve_data(OQP_td_mrsf_density, TA_TYPE_REAL64, nbf*nbf*7, (/7, nbf, nbf /), comment=OQP_td_mrsf_density)
    call infos%dat%reserve_data(OQP_td_p, TA_TYPE_REAL64, nbf_tri*2, (/ nbf_tri, 2 /), comment=OQP_td_p)
    call infos%dat%reserve_data(OQP_td_abxc, TA_TYPE_REAL64, nbf*nbf, (/ nbf, nbf /), comment=OQP_td_abxc)

    call data_has_tags(infos%dat, tags_alloc, module_name, subroutine_name, WITH_ABORT)
    call tagarray_get_data(infos%dat, OQP_WAO, wao)
    call tagarray_get_data(infos%dat, OQP_td_mrsf_density, td_mrsf_den)
    call tagarray_get_data(infos%dat, OQP_td_p, td_p)
    call tagarray_get_data(infos%dat, OQP_td_abxc, td_abxc)

    call data_has_tags(infos%dat, tags_required, module_name, subroutine_name, WITH_ABORT)
    call tagarray_get_data(infos%dat, OQP_FOCK_A, fock_a)
    call tagarray_get_data(infos%dat, OQP_FOCK_B, fock_b)
    call tagarray_get_data(infos%dat, OQP_E_MO_A, mo_energy_a)
    call tagarray_get_data(infos%dat, OQP_VEC_MO_A, mo_a)
    call tagarray_get_data(infos%dat, OQP_VEC_MO_B, mo_b)
    call tagarray_get_data(infos%dat, OQP_td_bvec_mo, bvec_mo)
    call tagarray_get_data(infos%dat, OQP_td_t, td_t)
    call tagarray_get_data(infos%dat, OQP_td_energies, mrsf_energies)

    ta => td_t(:,1)
    tb => td_t(:,2)

    target_state = min(infos%tddft%target_state, infos%tddft%nstate)
    if (target_state /=infos%tddft%target_state) then
      write(*,'(/1x,66("-")&
               &/1x,"WARNING: Target state has been changed to the max available nstates"/&
               &/1x,66("-")/)')
    end if

    ! Determine solver name for output
    if (infos%tddft%z_solver == 1) then
      solver_name = "GMRES"
    else
      solver_name = "CG"
    end if

    ! Save unrelaxed density matrices and the `b=A*x` vector for target state
    if (mrst==1 .or. mrst==3 ) then
      call mrsfxvec(infos, bvec_mo(:,target_state), bvec_mo_d(:,1))
      call sfdmat(bvec_mo_d(:,1), td_abxc, mo_a, ta, tb, nocca, noccb)
    else if (mrst==5 ) then
      call sfdmat(bvec_mo(:,target_state), td_abxc, mo_a, tb, ta, noccb, nocca)
    end if

    bvec(1:nbf,1:nbf,1:1) => td_abxc

  ! Initialize ERI calculations
    call int2_driver%init(basis, infos)
    call int2_driver%set_screening()

    write(*,'(/1x,71("-")&
             &/19x,"MRSF-DFT ENERGY GRADIENT CALCULATION"&
             &/1x,71("-")/)')

    write(iw,fmt='(5x,a/&
                  &5x,16("-")/&
                  &5x,a,x,i0,x,f17.10,x,"Hartree"/&
                  &5x,a,x,i0/&
                  &5x,a,x,e10.4/&
                  &5x,a,x,i0/&
                  &5x,a,x,a)') &
        'Z-vector options' &
      , 'Target state       is', target_state, infos%mol_energy%energy+mrsf_energies(target_state) &
      , 'Multiplicity       is', infos%tddft%mult &
      , 'Convergence        is', infos%tddft%zvconv &
      , 'Maximum iterations is', infos%control%maxit_zv &
      , 'Solver method      is', trim(solver_name)
    call flush(iw)

  ! Prepare for ROHF
    ! Fock matrices A and B
    if( roref )then
        wrk1t(1:nbf*nbf) => wrk1
  !   Alapha
      call orthogonal_transform_sym(nbf, nbf, fock_a, mo_a, nbf, wrk1)
      call unpack_matrix(wrk1t, fa)

  !   Beta
      call orthogonal_transform_sym(nbf, nbf, fock_b, mo_b, nbf, wrk1)
      call unpack_matrix(wrk1t, fb)
    end if

  ! Make density like part
    call unpack_matrix(ta, pa(:,:,1))
    call unpack_matrix(tb, pa(:,:,2))

  ! Initialize ERI calculations
    scale_exch = 1.0_dp
    scale_exch2 = 1.0_dp
    if (dft) then
       scale_exch = infos%dft%HFscale    !> Reference HF exchange
       scale_exch2 = infos%tddft%HFscale !> Response HF exchange
    end if

    if (mrst==1 .or. mrst==3 ) then

      int2_data_st = int2_mrsf_data_t( &
          d3 = fmrst1, &
          tamm_dancoff = .true., &
          scale_exchange = scale_exch2, &
          scale_coulomb = scale_exch2)

    else if( mrst==5  )then

      int2_data_q = int2_td_data_t( &
          d2=bvec, &
          int_apb = .false., &
          int_amb = .false., &
          tamm_dancoff = .true., &
          scale_exchange = scale_exch2)

    end if

    int2_data = int2_tdgrd_data_t( &
        d2 = pa, &
        int_apb = .true., &
        int_amb = .false., &
        tamm_dancoff = .false., &
        scale_exchange = scale_exch)

    call int2_driver%run(int2_data, &
            cam=dft.and.infos%dft%cam_flag, &
            alpha=infos%dft%cam_alpha, &
            beta=infos%dft%cam_beta,&
            mu=infos%dft%cam_mu)
    ab1 => int2_data%apb(:,:,:,1)

    pa = pa*2
    call utddft_fxc( &
        basis = basis, &
        molGrid = molGrid, &
        isVecs = .true., &
        wfa = mo_a, &
        wfb = mo_b, &
        fxa = ab1(:,:,1:1), &
        fxb = ab1(:,:,2:2), &
        dxa = pa(:,:,1:1), &
        dxb = pa(:,:,2:2), &
        nmtx = 1, &
        threshold = 1.0d-15, &
        infos = infos)

!   ALPHA: AO(M,N) -> MO(IA+)
    call mntoia(ab1(:,:,1), ab1_mo_a, mo_a, mo_a, nocca, nocca)

    call mntoia(ab1(:,:,2), ab1_mo_b, mo_b, mo_b, noccb, noccb)

    if (mrst==1 .or. mrst==3) then

      call iatogen(bvec_mo(:,target_state), wrk1, nocca, noccb)
      call mrsfcbc(infos, mo_a, mo_a, wrk1, fmrst1(1,:,:,:))

      fmrst1(1,7,:,:) = td_abxc

      td_mrsf_den(1:7,:,:) = fmrst1(1,1:7,:,:)

    ! Initialize ERI calculations
      call int2_driver%run(int2_data_st, &
            cam = dft.and.infos%dft%cam_flag, &
            alpha = infos%tddft%cam_alpha, &
            alpha_coulomb = infos%tddft%cam_alpha, &
            beta = infos%tddft%cam_beta,&
            beta_coulomb = infos%tddft%cam_beta, &
            mu = infos%tddft%cam_mu)
      fmrst2 => int2_data_st%f3(:,:,:,:,1)! ado2v, ado1v, adco1, adco2, ao21v, aco12, agdlr

    ! Scaling factor if triplet
      if (mrst==3) fmrst2(:,1:6,:,:) = -1.0_dp*fmrst2(:,1:6,:,:)

      ! Spin pair coupling
      if (infos%tddft%spc_coco /= infos%tddft%hfscale) &
         fmrst2(:,6,:,:) = fmrst2(:,6,:,:) * infos%tddft%spc_coco / infos%tddft%hfscale
      if (infos%tddft%spc_ovov /= infos%tddft%hfscale) &
         fmrst2(:,5,:,:) = fmrst2(:,5,:,:) * infos%tddft%spc_ovov / infos%tddft%hfscale
      if (infos%tddft%spc_coov /= infos%tddft%hfscale) &
         fmrst2(:,1:4,:,:) = fmrst2(:,1:4,:,:) * infos%tddft%spc_coov / infos%tddft%hfscale

      call orthogonal_transform('n', nbf, mo_a, fmrst2(1,7,:,:), wrk2, wrk1)

      call mrsfxvec(infos, bvec_mo(:,target_state), bvec_mo_d(:,1))

      call iatogen(bvec_mo_d(:,1), wrk3, nocca, noccb)

      call dgemm('n', 't', nbf, nocca, nbf, &
                 2.0_dp, wrk2, nbf, &
                         wrk3, nbf, &
                 0.0_dp, hxa, nbf)
      call dgemm('t', 'n', nbf, nbf, nocca, &
                 2.0_dp, wrk2, nbf, &
                         wrk3, nbf, &
                 0.0_dp, hxb, nbf)

   ! spin pair ov-ov, co-co, co-ov coupling
      call mrsfsp(hxa, hxb, mo_a, mo_a, wrk3, fmrst2(1,:,:,:), nocca, noccb)

   !  Unrelaxed difference density matries T_ij and T_ab
   !  Ta(i+,j+):= -X(i+,a-)*X(j+,a-) for singlet and triplet
      call dgemm('n', 't', nocca, nocca, nvirb, &
                -1.0_dp, bvec_mo_d, nocca, &
                         bvec_mo_d, nocca, &
                 0.0_dp, tij, nocca)

   !  Tb(a-,b-):= X(i+,a-)*X(i+,b-) for singlet and triplet
      call dgemm('t', 'n', nvirb, nvirb, nocca, &
                 1.0_dp, bvec_mo_d, nocca, &
                         bvec_mo_d, nocca, &
                 0.0_dp, tab, nvirb)

      call sfrorhs(rhs, hxa, hxb, ab1_mo_a, ab1_mo_b, &
                   Tij, Tab, Fa, Fb, nocca, noccb)

    else if(mrst==5) then

   !  Initialize ERI calculations
      call int2_driver%run(int2_data_q, &
            cam=dft.and.infos%dft%cam_flag, &
            alpha=infos%tddft%cam_alpha, &
            beta=infos%tddft%cam_beta,&
            mu=infos%tddft%cam_mu)

      call orthogonal_transform('n', nbf, mo_a, int2_data_q%amb(:,:,1,1), wrk2, wrk1)

      call iatogen(bvec_mo(:,target_state),wrk3,noccb,nocca)

      call dgemm('t', 'n', nbf, nbf, noccb, &
                 2.0_dp, wrk2, nbf, &
                         wrk3, nbf, &
                 0.0_dp, hxa, nbf)
      call dgemm('n', 't', nbf, noccb, nbf, &
                 2.0_dp, wrk2, nbf, &
                         wrk3, nbf, &
                 0.0_dp, hxb, nbf)

   !  Unrelaxed difference density matries T_ij and T_ab
   !  Ta(i+,j+):= -X(i+,a-)*X(j+,a-) for singlet and triplet
      call dgemm('n', 't', noccb, noccb, nvira, &
                -1.0_dp, bvec_mo(:,target_state), noccb, &
                         bvec_mo(:,target_state), noccb, &
                 0.0_dp, tij, noccb)

   !  Tb(a-,b-):= X(i+,a-)*X(i+,b-) for singlet and triplet
      call dgemm('t', 'n', nvira, nvira, noccb, &
                 1.0_dp, bvec_mo(:,target_state), noccb, &
                         bvec_mo(:,target_state), noccb, &
                 0.0_dp, tab, nvira)

      call mrsfqrorhs(rhs, hxa, hxb, ab1_mo_a, ab1_mo_b, &
                      tab, tij, fa, fb, nocca, noccb)
    end if

    write(*,'(/3x,25("-")&
             &/6x,"START Z-VECTOR LOOP (",A,")"&
             &/3x,25("-")/)') trim(solver_name)
    call flush(iw)

    call sfromcal(xm, xminv, mo_energy_a, fa, fb, nocca, noccb)

    ! Choose solver based on input option (0=CG, 1=GMRES)
    if (infos%tddft%z_solver == 1) then
      
      ! ============================================
      ! GMRES SOLVER
      ! ============================================
      
      ! Check preconditioner condition
      if (any(abs(xminv) < 1.0d-12) .or. any(abs(xminv) > 1.0d12)) then
        write(iw,'(" Warning: Preconditioner poorly conditioned, applying regularization")')
        where(abs(xminv) < 1.0d-12) xminv = sign(1.0d-12, xminv)
        where(abs(xminv) > 1.0d12) xminv = sign(1.0d12, xminv)
      end if
      
      ! Initial guess with same strategy as CG
      xk = 0.0_dp

      ! Call GMRES solver
      call gmres_solve( &
          apply_operator = apply_z_operator, &
          apply_precond = lambda_precond, &
          b = rhs, &
          x = xk, &
          n = lzdim, &
          restart = min(infos%tddft%gmres_dim, lzdim), &
          max_iter = infos%control%maxit_zv, &
          tol = cnvtol, &
          infos = infos, &
          basis = basis, &
          molGrid = molGrid, &
          int2_driver = int2_driver, &
          nocca = nocca, &
          noccb = noccb, &
          nbf = nbf, &
          mo_a = mo_a, &
          mo_b = mo_b, &
          mo_energy_a = mo_energy_a, &
          fa = fa, &
          fb = fb, &
          scale_exch = scale_exch, &
          dft = dft, &
          error_out = error, &
          iter_out = gmres_iter, &
          iw = iw)

      write(iw,'(/," Final Summary:")')
      write(iw,'(" GMRES total iterations: ", I4)') gmres_iter
      write(iw,'(" Final error norm      : ", 1p,e13.6)') error
      write(iw,'(" Convergence criterion : ", 1p,e13.6)') cnvtol
      call flush(iw)
      
      ! Clean up GMRES work arrays
      call cleanup_gmres_work()
      
    else
      
      ! ============================================
      ! ORIGINAL CONJUGATE GRADIENT SOLVER
      ! ============================================
      
      ! Initial guess with same strategy as CG
      xk = 0.0_dp

      call sfrogen(wrk1, wrk2, xk, nocca, noccb)
      ! Alpha
      call orthogonal_transform('t', nbf, mo_a, wrk1, pa(:,:,1), wrk3)
      ! Beta
      call orthogonal_transform('t', nbf, mo_b, wrk2, pa(:,:,2), wrk3)

      !****** INITIAL (A+B)*PK *************************************************
      ! Initialize ERI calculations
      call int2_data%clean()
      deallocate(int2_data)
      int2_data = int2_td_data_t( &
          d2 = pa, &
          int_apb = .true., &
          int_amb = .true., &
          tamm_dancoff = .false., &
          scale_exchange = scale_exch)

      call int2_driver%run(int2_data, &
              cam=dft.and.infos%dft%cam_flag, &
              alpha=infos%dft%cam_alpha, &
              beta=infos%dft%cam_beta,&
              mu=infos%dft%cam_mu)
      ab1 => int2_data%apb(:,:,:,1)
      ab2 => int2_data%amb(:,:,:,1)

      pa = pa*2
      call utddft_fxc( &
          basis = basis, &
          molGrid = molGrid, &
          isVecs = .true., &
          wfa = mo_a, &
          wfb = mo_b, &
          fxa = ab1(:,:,1:1), &
          fxb = ab1(:,:,2:2), &
          dxa = pa(:,:,1:1), &
          dxb = pa(:,:,2:2), &
          nmtx = 1, &
          threshold = 1.0d-15, &
          infos = infos)

      !   ALPHA: AO(M,N) -> MO(IA+) ... LPTMOA
      call mntoia(ab1(:,:,1), ab1_mo_a, mo_a, mo_a, nocca, nocca)

      wrk1 = 2*ab1(:,:,2) + ab2(:,:,2)
      call mntoia(wrk1, ab1_mo_b, mo_b, mo_b, noccb, noccb)

      call sfrolhs(lhs, xk, mo_energy_a, fa, fb, ab1_mo_a, ab1_mo_b, &
                   nocca, noccb)

      call pcgrbpini(errv, pk, error, rhs, xminv, lhs)

      write(iw,'(" Initial error =",3x,1p,e10.3,1x,"/",1p,e10.3)') error, cnvtol
      call flush(iw)

      ! -----------------------------------------------

      do iter = 1, infos%control%maxit_zv

        call sfrogen(wrk1, wrk2, pk, nocca, noccb)
        !     Alpha
        call orthogonal_transform('t', nbf, mo_a, wrk1, pa(:,:,1), wrk3)
        !     Beta
        call orthogonal_transform('t', nbf, mo_b, wrk2, pa(:,:,2), wrk3)

        !     (A+B)*PK
        call int2_data%clean()
        deallocate(int2_data)
        int2_data = int2_tdgrd_data_t( &
            d2 = pa, &
            int_apb = .true., &
            int_amb = .false., &
            tamm_dancoff = .false., &
            scale_exchange = scale_exch)

        call int2_driver%run(int2_data, &
              cam=dft.and.infos%dft%cam_flag, &
              alpha=infos%dft%cam_alpha, &
              beta=infos%dft%cam_beta,&
              mu=infos%dft%cam_mu)
        ab1 => int2_data%apb(:,:,:,1)

        !ab1 = ab1/2
        call symmetrize_matrix(pa(:,:,1), nbf)
        call symmetrize_matrix(pa(:,:,2), nbf)
        call utddft_fxc( &
            basis = basis, &
            molGrid = molGrid, &
            isVecs = .true., &
            wfa = mo_a, &
            wfb = mo_b, &
            fxa = ab1(:,:,1:1), &
            fxb = ab1(:,:,2:2), &
            dxa = pa(:,:,1:1), &
            dxb = pa(:,:,2:2), &
            nmtx = 1, &
            threshold = 1.0d-15, &
            infos = infos)

        !     ALPHA: AO(M,N) -> MO(IA+) ... LPTMOA
        call mntoia(ab1(:,:,1), ab1_mo_a, mo_a, mo_a, nocca, nocca)

        call mntoia(ab1(:,:,2), ab1_mo_b, mo_b, mo_b, noccb, noccb)

        call sfrolhs(lhs, pk, mo_energy_a, fa, fb, ab1_mo_a, ab1_mo_b, &
                     nocca, noccb)

        alpha = 1.0_dp/dot_product(pk, lhs)

        xk = xk + pk * alpha
        errv = errv - alpha*lhs

        error = dot_product(errv, errv)
        write(iw,'(" Iter#",I2," Error =",&
              &3x,1p,e10.3,1x,"/",1p,e10.3)') &
                iter, error, cnvtol
        call flush(iw)

        if (error<cnvtol) exit

        call pcgb(pk, errv, xminv)

      end do
      
    end if  ! End solver selection

! -----------------------------------------------
    if (error>cnvtol) then
       infos%mol_energy%Z_Vector_converged=.false.
       write(*,'(/3x,24("-")&
             &/6x,"Z-Vector not converged"&
             &/3x,24("-")/)')
    else
       infos%mol_energy%Z_Vector_converged=.true.
       write(*,'(/3x,24("-")&
             &/6x,"Z-Vector converged"&
             &/3x,24("-")/)')
    endif

    call flush(iw)

    if (mrst==1 .or. mrst==3) then

      call sfropcal(wrk1, wrk2, tij, tab, xk, nocca, noccb)

    else if (mrst==5) then

      call mrsfqropcal(wrk1, wrk2, tab, tij, xk, nocca, noccb)

    end if

 !  Update density for alpha
    call orthogonal_transform('t', nbf, mo_a, wrk1, pa(:,:,1), wrk3)

 !  Update density for beta
    call orthogonal_transform('t', nbf, mo_b, wrk2, pa(:,:,2), wrk3)
    call int2_data%clean()
    deallocate(int2_data)
    int2_data = int2_tdgrd_data_t( &
        d2 = pa, &
        int_apb = .true., &
        int_amb = .false., &
        tamm_dancoff = .false., &
        scale_exchange = scale_exch)

    call int2_driver%run(int2_data, &
            cam=dft.and.infos%dft%cam_flag, &
            alpha=infos%dft%cam_alpha, &
            beta=infos%dft%cam_beta,&
            mu=infos%dft%cam_mu)
    ab1 => int2_data%apb(:,:,:,1)

    call symmetrize_matrix(pa(:,:,1), nbf)
    call symmetrize_matrix(pa(:,:,2), nbf)
    call pack_matrix(pa(:,:,1), td_p(:,1))
    call pack_matrix(pa(:,:,2), td_p(:,2))

    td_p = 0.5_dp*td_p

    call utddft_fxc( &
        basis = basis, &
        molGrid = molGrid, &
        isVecs = .true., &
        wfa = mo_a, &
        wfb = mo_b, &
        fxa = ab1(:,:,1:1), &
        fxb = ab1(:,:,2:2), &
        dxa = pa(:,:,1:1), &
        dxb = pa(:,:,2:2), &
        nmtx = 1, &
        threshold = 1.0d-15, &
        infos = infos)

!   ALPHA AO(M,N) -> MO(I-,J-) ... LPPIJA
    call dgemm('n', 'n', nbf, nocca, nbf,  &
               1.0_dp, ab1(:,:,1), nbf,  &
                       mo_a, nbf,  &
               0.0_dp, wrk2, nbf)
    call dgemm('t', 'n', nocca, nocca, nbf,  &
               1.0_dp, mo_a,  nbf,  &
                       wrk2,  nbf,  &
               0.0_dp, ppija, nocca)
!   BETA: AO(M,N) -> MO(I-,J-) ... LPPIJB
    call dgemm('n', 'n', nbf, noccb, nbf,  &
               1.0_dp, ab1(:,:,2), nbf,  &
                       mo_b, nbf,  &
               0.0_dp, wrk2, nbf)
    call dgemm('t', 'n', noccb, noccb, nbf,  &
               1.0_dp, mo_b,  nbf,  &
                       wrk2,  nbf,  &
               0.0_dp, ppijb, noccb)

!   Calculate W (in MO basis)
    wmo => wrk3
    wmo = 0
    if (mrst==1 .or. mrst==3) then

      call mrsfrowcal(wmo, mo_energy_a, fa, fb, xk, &
                      hxa, hxb, ppija, ppijb, &
                      nocca, noccb)

    else if (mrst==5) then

      call mrsfqrowcal(wmo, mo_energy_a, fa, fb, xk, &
                       hxa, hxb, ppija, ppijb, &
                       nocca, noccb)

    end if

    call orthogonal_transform('t', nbf, mo_a, wmo, wrk2, wrk1)
    call symmetrize_matrix(wrk2, nbf)
    call pack_matrix(wrk2, wao)
    wao = wao*0.5_dp
!   ROHF, half one more time:
    wao = wao*0.5_dp

    call int2_driver%clean()

    if (dft) call dftclean(infos)

    call measure_time(print_total=1, log_unit=iw)
    ! Clean up GMRES work arrays
    call cleanup_gmres_work()

    close(iw)

  contains

    ! Lambda wrapper for preconditioner
    subroutine lambda_precond(x_in, x_out)
      real(kind=dp), intent(in) :: x_in(:)
      real(kind=dp), intent(out) :: x_out(:)
      call apply_z_precond(x_in, x_out, xminv)
    end subroutine lambda_precond

  end subroutine tdhf_mrsf_z_vector

end module tdhf_mrsf_z_vector_mod