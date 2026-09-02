module tdhf_hessian_response_mod

  use precision, only: dp
  use, intrinsic :: iso_fortran_env, only: error_unit
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

  implicit none

  private
  public :: solve_tdhf_amplitude_response
  public :: solve_tdhf_z_response
  public :: solve_mrsf_z_response_matrix_free
  public :: solve_mrsf_z_response_batch_matrix_free
  public :: complete_rhf_orbital_response
  public :: assemble_tdhf_sigma_derivative
  public :: solve_mrsf_tda_response_dense
  public :: solve_mrsf_tda_response_matrix_free
  public :: solve_mrsf_tda_response_batch_matrix_free
  public :: solve_mrsf_tda_cluster_response_matrix_free
  public :: assemble_mrsf_tda_eigenvalue_hessian

  abstract interface
    subroutine mrsf_tda_operator(vector, result, status)
      import dp
      real(kind=dp), intent(in) :: vector(:)
      real(kind=dp), intent(out) :: result(:)
      integer, intent(out) :: status
    end subroutine mrsf_tda_operator

    subroutine mrsf_tda_batch_operator(vectors, results, status)
      import dp
      real(kind=dp), intent(in) :: vectors(:,:)
      real(kind=dp), intent(out) :: results(:,:)
      integer, intent(out) :: status
    end subroutine mrsf_tda_batch_operator
  end interface

contains

!###############################################################################

  subroutine assemble_mrsf_tda_eigenvalue_hessian(x0,dax,dx,d2a_expect, &
                                                    hessian,asymmetry,status)
    ! Assemble the electronic-state eigenvalue part of the MRSF-TDA Hessian,
    !
    ! d2omega(K,L) = X^T d2A(K,L) X + 2 (dX(L))^T dA(K) X.
    !
    ! dax(:,K) is the differentiated response-operator action dA(K)X and
    ! d2a_expect(K,L) is the explicit fixed-X expectation X^T d2A(K,L)X.
    ! For an exact response the raw matrix is symmetric.  We report its
    ! antisymmetric diagnostic before returning the symmetrized Hessian.

    real(kind=dp), intent(in) :: x0(:),dax(:,:),dx(:,:),d2a_expect(:,:)
    real(kind=dp), intent(out) :: hessian(:,:),asymmetry
    integer, intent(out) :: status

    real(kind=dp), allocatable :: raw(:,:)
    integer :: k,l,n,ncoord

    n = size(x0)
    ncoord = size(dax,2)
    status = 0
    hessian = 0.0_dp
    asymmetry = 0.0_dp
    if (n <= 0 .or. ncoord <= 0 .or. size(dax,1) /= n .or. &
        any(shape(dx) /= [n,ncoord]) .or. &
        any(shape(d2a_expect) /= [ncoord,ncoord]) .or. &
        any(shape(hessian) /= [ncoord,ncoord])) then
      status = -1
      return
    end if
    if (abs(dot_product(x0,x0)-1.0_dp) > 1.0e-10_dp .or. &
        maxval(abs(matmul(x0,dx))) > 1.0e-10_dp) then
      status = -2
      return
    end if
    allocate(raw(ncoord,ncoord))
    do l=1,ncoord
      do k=1,ncoord
        raw(k,l) = d2a_expect(k,l)+2.0_dp*dot_product(dx(:,l),dax(:,k))
      end do
    end do
    asymmetry = maxval(abs(raw-transpose(raw)))
    hessian = 0.5_dp*(raw+transpose(raw))
    deallocate(raw)
  end subroutine assemble_mrsf_tda_eigenvalue_hessian

!###############################################################################

  subroutine solve_mrsf_tda_response_dense(a, omega, x0, dax, dx, domega, &
                                             residual_max, status, pivot_tol)
    ! Dense reference solution of the differentiated symmetric MRSF-TDA
    ! eigenproblem for an isolated state,
    !
    !   (A-omega I) dX = -(dA-domega I) X,   X^T dX = 0.
    !
    ! The bordered equation fixes the otherwise undetermined component along
    ! X.  It is intentionally retained as an independent small-system
    ! reference.  Molecular calculations use the same projected equation with
    ! a matrix-vector product and an indefinite iterative solver.

    real(kind=dp), intent(in) :: a(:,:), omega, x0(:), dax(:,:)
    real(kind=dp), intent(out) :: dx(:,:), domega(:), residual_max
    integer, intent(out) :: status
    real(kind=dp), intent(in), optional :: pivot_tol

    real(kind=dp), allocatable :: bordered(:,:), rhs(:), solution(:), check(:)
    real(kind=dp) :: norm2, threshold
    integer :: coord, n, ncoord, solve_status

    n = size(x0)
    ncoord = size(dax,2)
    status = 0
    residual_max = 0.0_dp
    dx = 0.0_dp
    domega = 0.0_dp
    threshold = 1.0e-12_dp
    if (present(pivot_tol)) threshold = pivot_tol

    if (n <= 0 .or. any(shape(a) /= [n,n]) .or. size(dax,1) /= n .or. &
        any(shape(dx) /= [n,ncoord]) .or. size(domega) /= ncoord .or. &
        threshold <= 0.0_dp) then
      status = -1
      return
    end if
    norm2 = dot_product(x0,x0)
    if (abs(norm2-1.0_dp) > 1.0e-10_dp .or. &
        maxval(abs(matmul(a,x0)-omega*x0)) > 1.0e-9_dp .or. &
        maxval(abs(a-transpose(a))) > 1.0e-10_dp) then
      status = -2
      return
    end if

    allocate(bordered(n+1,n+1),rhs(n+1),solution(n+1),check(n))
    bordered = 0.0_dp
    bordered(1:n,1:n) = a
    do coord = 1,n
      bordered(coord,coord) = bordered(coord,coord)-omega
    end do
    bordered(1:n,n+1) = x0
    bordered(n+1,1:n) = x0

    do coord = 1,ncoord
      domega(coord) = dot_product(x0,dax(:,coord))
      rhs(1:n) = -dax(:,coord)+domega(coord)*x0
      rhs(n+1) = 0.0_dp
      call solve_linear_pivot(bordered,rhs,solution,threshold,solve_status)
      if (solve_status /= 0) then
        status = coord
        exit
      end if
      dx(:,coord) = solution(1:n)
      check = matmul(a,dx(:,coord))-omega*dx(:,coord) &
        +dax(:,coord)-domega(coord)*x0
      residual_max = max(residual_max,maxval(abs(check)), &
        abs(dot_product(x0,dx(:,coord))))
    end do
    deallocate(bordered,rhs,solution,check)

  end subroutine solve_mrsf_tda_response_dense

!###############################################################################

  subroutine solve_mrsf_tda_response_matrix_free(apply_operator,omega,x0,dax, &
      dx,domega,residual_max,status,tol,maxit,restart)
    ! Matrix-free differentiated MRSF-TDA response for an isolated state.
    ! The physical response equation is solved only in the complement of X,
    !
    !   P(A-omega I)P dX = -P(dA)X,       P=I-XX^T,
    !
    ! and the X direction is given a unit eigenvalue to make the full-space
    ! Krylov equation nonsingular.  Restarted GMRES is used because the
    ! projected operator is symmetric but generally indefinite for any root
    ! above the lowest state.  No explicit response matrix is constructed.

    procedure(mrsf_tda_operator) :: apply_operator
    real(kind=dp), intent(in) :: omega,x0(:),dax(:,:)
    real(kind=dp), intent(out) :: dx(:,:),domega(:),residual_max
    integer, intent(out) :: status
    real(kind=dp), intent(in), optional :: tol
    integer, intent(in), optional :: maxit,restart

    real(kind=dp), allocatable :: rhs(:),solution(:),check(:),ax(:)
    real(kind=dp) :: norm2,solve_tol
    integer :: coordinate,n,ncoord,niter,nrestart,operator_status,solve_status

    n=size(x0); ncoord=size(dax,2)
    status=0; residual_max=0.0_dp; dx=0.0_dp; domega=0.0_dp
    solve_tol=1.0e-11_dp
    if(present(tol)) solve_tol=tol
    niter=max(100,4*n)
    if(present(maxit)) niter=maxit
    nrestart=min(max(8,n/4),max(1,n))
    if(present(restart)) nrestart=restart
    if(n<=0 .or. ncoord<=0 .or. size(dax,1)/=n .or. &
       any(shape(dx)/=[n,ncoord]) .or. size(domega)/=ncoord .or. &
       solve_tol<=0.0_dp .or. niter<=0 .or. nrestart<=0) then
      status=-1
      return
    end if
    norm2=dot_product(x0,x0)
    if(abs(norm2-1.0_dp)>1.0e-10_dp) then
      status=-2
      return
    end if
    allocate(rhs(n),solution(n),check(n),ax(n))
    call apply_operator(x0,ax,operator_status)
    residual_max=maxval(abs(ax-omega*x0))
    if(operator_status/=0 .or. residual_max>1.0e-8_dp) then
      write(error_unit,'(A,I0,A,1P,E15.7,A,E15.7)') &
        ' MRSF response eigenpair check: operator status ',operator_status, &
        ', max residual = ',residual_max,', omega = ',omega
      status=-3
      deallocate(rhs,solution,check,ax)
      return
    end if
    do coordinate=1,ncoord
      domega(coordinate)=dot_product(x0,dax(:,coordinate))
      rhs=-dax(:,coordinate)+domega(coordinate)*x0
      rhs=rhs-x0*dot_product(x0,rhs)
      solution=0.0_dp
      call solve_projected_gmres(apply_operator,omega,x0,rhs,solution, &
        solve_tol,niter,nrestart,solve_status)
      if(solve_status/=0) then
        write(error_unit,'(A,I0,A,I0)') &
          'MRSF projected amplitude solve failed for coordinate ', &
          coordinate,' with status ',solve_status
        status=coordinate
        exit
      end if
      solution=solution-x0*dot_product(x0,solution)
      dx(:,coordinate)=solution
      call apply_operator(solution,ax,operator_status)
      if(operator_status/=0) then
        status=coordinate
        exit
      end if
      check=ax-omega*solution+dax(:,coordinate)-domega(coordinate)*x0
      residual_max=max(residual_max,maxval(abs(check)), &
        abs(dot_product(x0,solution)))
    end do
    deallocate(rhs,solution,check,ax)
  end subroutine solve_mrsf_tda_response_matrix_free

!###############################################################################

  subroutine solve_mrsf_tda_response_batch_matrix_free(apply_operator,omega, &
      x0,dax,dx,domega,residual_max,status,tol,maxit,apply_preconditioner)
    ! Matrix-free multi-right-hand-side response in the physical MRSF-TDA
    ! space.  This follows the efficient organization of the GAMESS analytic
    ! TDHF/TDDFT Hessian: all nuclear perturbations share each expensive AO
    ! response pass.  All right-hand sides contribute trial vectors to one
    ! shared, restarted residual-minimization subspace.  Thus related nuclear
    ! perturbations reuse both the AO response pass and each other's search
    ! directions.  No response matrix or determinant expansion is formed.

    procedure(mrsf_tda_batch_operator) :: apply_operator
    procedure(mrsf_tda_batch_operator), optional :: apply_preconditioner
    real(kind=dp), intent(in) :: omega,x0(:),dax(:,:)
    real(kind=dp), intent(out) :: dx(:,:),domega(:),residual_max
    integer, intent(out) :: status
    real(kind=dp), intent(in), optional :: tol
    integer, intent(in), optional :: maxit

    real(kind=dp), allocatable :: rhs(:,:),trial(:,:),applied(:,:),check(:,:), &
      projected(:,:),operator_image(:,:),norm_rhs(:),parallel(:), &
      parallel_result(:)
    real(kind=dp) :: solve_tol,subspace_residual
    integer :: coordinate,n,ncoord,niter,operator_status,solve_status

    n=size(x0)
    ncoord=size(dax,2)
    status=0
    residual_max=0.0_dp
    dx=0.0_dp
    domega=0.0_dp
    solve_tol=1.0e-10_dp
    if(present(tol)) solve_tol=tol
    niter=max(200,4*n)
    if(present(maxit)) niter=maxit
    if(n<=0 .or. ncoord<=0 .or. size(dax,1)/=n .or. &
       any(shape(dx)/=[n,ncoord]) .or. size(domega)/=ncoord .or. &
       .not.ieee_is_finite(omega) .or. &
       .not.ieee_is_finite(solve_tol) .or. solve_tol<=0.0_dp .or. &
       niter<=0 .or. any(.not.ieee_is_finite(x0)) .or. &
       any(.not.ieee_is_finite(dax))) then
      status=-1
      return
    end if
    if(abs(dot_product(x0,x0)-1.0_dp)>1.0e-10_dp) then
      status=-2
      return
    end if

    allocate(rhs(n,ncoord),trial(n,ncoord),applied(n,ncoord), &
      check(n,ncoord),projected(n,ncoord),operator_image(n,ncoord), &
      norm_rhs(ncoord),parallel(ncoord),parallel_result(ncoord))

    ! Validate the supplied eigenpair through the same production batch
    ! operator that will be used by every response direction.
    trial=0.0_dp
    trial(:,1)=x0
    call apply_operator(trial(:,1:1),applied(:,1:1),operator_status)
    if(operator_status/=0 .or. any(.not.ieee_is_finite(applied(:,1)))) then
      status=-3
      go to 900
    end if
    residual_max=maxval(abs(applied(:,1)-omega*x0))
    if(residual_max>1.0e-8_dp) then
      write(error_unit,'(A,1P,E15.7,A,E15.7)') &
        ' MRSF batch response eigenpair max residual = ',residual_max, &
        ', omega = ',omega
      status=-3
      go to 900
    end if

    do coordinate=1,ncoord
      domega(coordinate)=dot_product(x0,dax(:,coordinate))
      rhs(:,coordinate)=-dax(:,coordinate)+domega(coordinate)*x0
      rhs(:,coordinate)=rhs(:,coordinate)- &
        x0*dot_product(x0,rhs(:,coordinate))
      norm_rhs(coordinate)=sqrt(dot_product(rhs(:,coordinate), &
        rhs(:,coordinate)))
    end do
    if(any(.not.ieee_is_finite(norm_rhs))) then
      status=-4
      go to 900
    end if

    if(present(apply_preconditioner)) then
      call solve_shared_block_subspace(apply_projected_batch,rhs,dx, &
        subspace_residual,solve_status,solve_tol,niter, &
        apply_preconditioner=apply_projected_preconditioner)
    else
      call solve_shared_block_subspace(apply_projected_batch,rhs,dx, &
        subspace_residual,solve_status,solve_tol,niter)
    end if
    if(solve_status/=0) then
      status=solve_status
      go to 900
    end if

    ! One final batched application certifies every returned response against
    ! the unprojected differentiated eigenvalue equation.
    call apply_operator(dx,applied,operator_status)
    if(operator_status/=0 .or. any(.not.ieee_is_finite(applied))) then
      status=-5
      go to 900
    end if
    do coordinate=1,ncoord
      check(:,coordinate)=applied(:,coordinate)-omega*dx(:,coordinate)+ &
        dax(:,coordinate)-domega(coordinate)*x0
      residual_max=max(residual_max,maxval(abs(check(:,coordinate))), &
        abs(dot_product(x0,dx(:,coordinate))))
      if(sqrt(dot_product(check(:,coordinate),check(:,coordinate)))> &
         max(1.0e-8_dp,100.0_dp*solve_tol)* &
           max(1.0_dp,norm_rhs(coordinate))) then
        status=coordinate
        go to 900
      end if
    end do

900 continue
    deallocate(rhs,trial,applied,check,projected,operator_image,norm_rhs, &
      parallel,parallel_result)

  contains

    subroutine apply_projected_batch(vectors,results,local_status)
      real(kind=dp), intent(in) :: vectors(:,:)
      real(kind=dp), intent(out) :: results(:,:)
      integer, intent(out) :: local_status
      integer :: nvec

      local_status=0
      results=0.0_dp
      nvec=size(vectors,2)
      if(size(vectors,1)/=n .or. any(shape(results)/=[n,nvec]) .or. &
         nvec<=0 .or. any(.not.ieee_is_finite(vectors))) then
        local_status=-1
        return
      end if
      parallel(1:nvec)=matmul(x0,vectors)
      projected(:,1:nvec)=vectors- &
        spread(x0,2,nvec)*spread(parallel(1:nvec),1,n)
      call apply_operator(projected(:,1:nvec),operator_image(:,1:nvec), &
        local_status)
      if(local_status/=0 .or. &
         any(.not.ieee_is_finite(operator_image(:,1:nvec)))) then
        local_status=-1
        results=0.0_dp
        return
      end if
      results=operator_image(:,1:nvec)-omega*projected(:,1:nvec)
      parallel_result(1:nvec)=matmul(x0,results)
      results=results- &
        spread(x0,2,nvec)*spread(parallel_result(1:nvec),1,n)+ &
        spread(x0,2,nvec)*spread(parallel(1:nvec),1,n)
      if(any(.not.ieee_is_finite(results))) local_status=-1
    end subroutine apply_projected_batch

    subroutine apply_projected_preconditioner(vectors,results,local_status)
      real(kind=dp), intent(in) :: vectors(:,:)
      real(kind=dp), intent(out) :: results(:,:)
      integer, intent(out) :: local_status
      integer :: nvec

      nvec=size(vectors,2)
      local_status=0
      if(size(vectors,1)/=n .or. any(shape(results)/=[n,nvec]) .or. &
         nvec<=0 .or. any(.not.ieee_is_finite(vectors))) then
        results=0.0_dp
        local_status=-1
        return
      end if
      call apply_preconditioner(vectors,results,local_status)
      if(local_status/=0 .or. any(.not.ieee_is_finite(results))) then
        results=0.0_dp
        local_status=-1
        return
      end if
      parallel_result(1:nvec)=matmul(x0,results)
      results=results-spread(x0,2,nvec)* &
        spread(parallel_result(1:nvec),1,n)
    end subroutine apply_projected_preconditioner

  end subroutine solve_mrsf_tda_response_batch_matrix_free

!###############################################################################

  subroutine solve_mrsf_tda_cluster_response_matrix_free(apply_operator, &
      cluster_energies,cluster_vectors,dax,response,effective_derivative, &
      residual_max,status,tol,maxit,restart)
    ! Gauge-fixed response of a near-degenerate MRSF-TDA root cluster.
    !
    ! Methodological starting point: Hiroya Nakata's analytical TDHF/TDDFT
    ! Hessian response formulation.  Here its isolated-root response equation
    ! is extended to an invariant cluster projector.  For P=X X^T, Q=I-P,
    ! and Hc=X^T A X, the returned complement response Z satisfies
    !
    !   Q (A Z-Z Hc) = -Q (dA X),        X^T Z=0,
    !   Weff = sym[X^T (dA X)].
    !
    ! Hc is evaluated from the operator callback rather than assumed diagonal.
    ! Consequently the equation remains covariant under X -> X U even when U
    ! rotates unequal but clustered roots: Hc -> U^T Hc U and Z -> Z U.
    ! cluster_energies supplies the invariant cluster spectrum diagnostic.
    !
    ! The callback is assumed to act on a complete, already-defined physical
    ! MRSF-TDA amplitude vector spanning the CO, OV, CV, and OO sectors.  Its
    ! production action may be assembled through seven AO density blocks; the
    ! density channels are not identified with the four amplitude sectors.
    ! This solver neither constructs nor interprets determinants or Slater
    ! strings; it acts only on the supplied spin-adapted physical vector space.

    procedure(mrsf_tda_operator) :: apply_operator
    real(kind=dp), intent(in) :: cluster_energies(:),cluster_vectors(:,:),dax(:,:)
    real(kind=dp), intent(out) :: response(:,:),effective_derivative(:,:),residual_max
    integer, intent(out) :: status
    real(kind=dp), intent(in), optional :: tol
    integer, intent(in), optional :: maxit,restart

    real(kind=dp), allocatable :: ax(:,:),cluster_operator(:,:),overlap(:,:), &
      projector_work(:,:),rhs_matrix(:,:),rhs(:),solution(:),check(:,:), &
      trial(:,:),trial_q(:,:),atrial(:,:),result_matrix(:,:),parallel(:,:)
    real(kind=dp), allocatable :: computed_energies(:),sorted_energies(:)
    real(kind=dp) :: energy_scale,operator_scale,solve_tol,spectrum_error
    integer :: column,eigensolver_status,n,ncluster,niter,nrestart, &
      operator_status,solve_status

    n=size(cluster_vectors,1)
    ncluster=size(cluster_vectors,2)
    status=0
    residual_max=0.0_dp
    response=0.0_dp
    effective_derivative=0.0_dp
    solve_tol=1.0e-11_dp
    if(present(tol)) solve_tol=tol
    niter=max(200,4*n*ncluster)
    if(present(maxit)) niter=maxit
    nrestart=min(max(12,n*ncluster/4),max(1,n*ncluster))
    if(present(restart)) nrestart=restart
    if(n<=0 .or. ncluster<=0 .or. ncluster>=n .or. &
       size(cluster_energies)/=ncluster .or. &
       any(shape(dax)/=[n,ncluster]) .or. &
       any(shape(response)/=[n,ncluster]) .or. &
       any(shape(effective_derivative)/=[ncluster,ncluster]) .or. &
       solve_tol<=0.0_dp .or. niter<=0 .or. nrestart<=0) then
      status=-1
      return
    end if

    allocate(ax(n,ncluster),cluster_operator(ncluster,ncluster), &
      overlap(ncluster,ncluster),projector_work(n,ncluster), &
      rhs_matrix(n,ncluster),rhs(n*ncluster),solution(n*ncluster), &
      check(n,ncluster),trial(n,ncluster),trial_q(n,ncluster), &
      atrial(n,ncluster),result_matrix(n,ncluster), &
      parallel(ncluster,ncluster),computed_energies(ncluster), &
      sorted_energies(ncluster))

    overlap=matmul(transpose(cluster_vectors),cluster_vectors)
    do column=1,ncluster
      overlap(column,column)=overlap(column,column)-1.0_dp
    end do
    if(maxval(abs(overlap))>1.0e-10_dp) then
      status=-2
      go to 900
    end if

    do column=1,ncluster
      call apply_operator(cluster_vectors(:,column),ax(:,column),operator_status)
      if(operator_status/=0) then
        status=-3
        go to 900
      end if
    end do
    cluster_operator=matmul(transpose(cluster_vectors),ax)
    cluster_operator=0.5_dp*(cluster_operator+transpose(cluster_operator))
    projector_work=ax-matmul(cluster_vectors,cluster_operator)
    operator_scale=max(1.0_dp,maxval(abs(ax)),maxval(abs(cluster_operator)))
    if(maxval(abs(projector_work))>max(1.0e-9_dp,100.0_dp*solve_tol)*operator_scale) then
      status=-4
      go to 900
    end if

    ! Compare spectra, rather than diagonal elements, so the diagnostic does
    ! not select a privileged basis within the root cluster.
    energy_scale=max(1.0_dp,maxval(abs(cluster_energies)))
    call symmetric_eigenvalues_jacobi(cluster_operator,computed_energies, &
      eigensolver_status)
    if(eigensolver_status/=0) then
      status=-5
      go to 900
    end if
    sorted_energies=cluster_energies
    call sort_real_ascending(sorted_energies)
    spectrum_error=maxval(abs(computed_energies-sorted_energies))
    if(spectrum_error>max(1.0e-8_dp,1000.0_dp*solve_tol)*energy_scale) then
      status=-5
      go to 900
    end if

    effective_derivative=matmul(transpose(cluster_vectors),dax)
    effective_derivative=0.5_dp*(effective_derivative+ &
      transpose(effective_derivative))
    rhs_matrix=-dax+matmul(cluster_vectors, &
      matmul(transpose(cluster_vectors),dax))
    rhs=reshape(rhs_matrix,[n*ncluster])
    solution=0.0_dp
    call solve_general_gmres(apply_cluster_projected_operator,rhs,solution, &
      solve_tol,niter,nrestart,solve_status)
    if(solve_status/=0) then
      status=10+abs(solve_status)
      go to 900
    end if
    response=reshape(solution,[n,ncluster])
    response=response-matmul(cluster_vectors, &
      matmul(transpose(cluster_vectors),response))

    do column=1,ncluster
      call apply_operator(response(:,column),atrial(:,column),operator_status)
      if(operator_status/=0) then
        status=-3
        go to 900
      end if
    end do
    check=atrial-matmul(response,cluster_operator)+dax- &
      matmul(cluster_vectors,effective_derivative)
    residual_max=max(maxval(abs(check)), &
      maxval(abs(matmul(transpose(cluster_vectors),response))))
    if(residual_max>max(1.0e-8_dp,100.0_dp*solve_tol)* &
       max(1.0_dp,maxval(abs(dax)))) status=-6

900 continue
    deallocate(ax,cluster_operator,overlap,projector_work,rhs_matrix,rhs, &
      solution,check,trial,trial_q,atrial,result_matrix,parallel, &
      computed_energies,sorted_energies)

  contains

    subroutine apply_cluster_projected_operator(vector,result,local_status)
      real(kind=dp), intent(in) :: vector(:)
      real(kind=dp), intent(out) :: result(:)
      integer, intent(out) :: local_status
      integer :: local_column

      local_status=0
      if(size(vector)/=n*ncluster .or. size(result)/=n*ncluster) then
        local_status=-1
        result=0.0_dp
        return
      end if
      trial=reshape(vector,[n,ncluster])
      parallel=matmul(transpose(cluster_vectors),trial)
      trial_q=trial-matmul(cluster_vectors,parallel)
      do local_column=1,ncluster
        call apply_operator(trial_q(:,local_column), &
          atrial(:,local_column),local_status)
        if(local_status/=0) then
          result=0.0_dp
          return
        end if
      end do
      result_matrix=atrial-matmul(trial_q,cluster_operator)
      result_matrix=result_matrix-matmul(cluster_vectors, &
        matmul(transpose(cluster_vectors),result_matrix)) &
        +matmul(cluster_vectors,parallel)
      result=reshape(result_matrix,[n*ncluster])
    end subroutine apply_cluster_projected_operator

  end subroutine solve_mrsf_tda_cluster_response_matrix_free

!###############################################################################

  subroutine symmetric_eigenvalues_jacobi(matrix,eigenvalues,status)
    ! Small, dependency-free eigensolver used only to validate the supplied
    ! cluster spectrum.  The response itself remains matrix free.
    real(kind=dp), intent(in) :: matrix(:,:)
    real(kind=dp), intent(out) :: eigenvalues(:)
    integer, intent(out) :: status

    real(kind=dp), allocatable :: work(:,:)
    real(kind=dp) :: akp,akq,app,apq,aqq,c,max_off,scale,s,t,tau
    integer :: i,iteration,j,k,max_iterations,n,p,q

    n=size(eigenvalues)
    status=0
    eigenvalues=0.0_dp
    if(n<=0 .or. any(shape(matrix)/=[n,n])) then
      status=-1
      return
    end if
    allocate(work(n,n))
    work=0.5_dp*(matrix+transpose(matrix))
    scale=max(1.0_dp,maxval(abs(work)))
    max_iterations=max(30,20*n*n)
    max_off=0.0_dp
    do iteration=1,max_iterations
      max_off=0.0_dp; p=1; q=1
      do j=2,n
        do i=1,j-1
          if(abs(work(i,j))>max_off) then
            max_off=abs(work(i,j)); p=i; q=j
          end if
        end do
      end do
      if(max_off<=100.0_dp*epsilon(1.0_dp)*scale) exit
      app=work(p,p); apq=work(p,q); aqq=work(q,q)
      tau=(aqq-app)/(2.0_dp*apq)
      t=sign(1.0_dp,tau)/(abs(tau)+hypot(1.0_dp,tau))
      c=1.0_dp/sqrt(1.0_dp+t*t)
      s=t*c
      do k=1,n
        if(k==p .or. k==q) cycle
        akp=work(k,p); akq=work(k,q)
        work(k,p)=c*akp-s*akq; work(p,k)=work(k,p)
        work(k,q)=s*akp+c*akq; work(q,k)=work(k,q)
      end do
      work(p,p)=app-t*apq
      work(q,q)=aqq+t*apq
      work(p,q)=0.0_dp; work(q,p)=0.0_dp
    end do
    if(max_off>100.0_dp*epsilon(1.0_dp)*scale) then
      status=-2
    else
      do i=1,n
        eigenvalues(i)=work(i,i)
      end do
      call sort_real_ascending(eigenvalues)
    end if
    deallocate(work)
  end subroutine symmetric_eigenvalues_jacobi

!###############################################################################

  subroutine sort_real_ascending(values)
    real(kind=dp), intent(inout) :: values(:)
    real(kind=dp) :: value
    integer :: i,j
    do i=2,size(values)
      value=values(i)
      j=i-1
      do while(j>=1)
        if(values(j)<=value) exit
        values(j+1)=values(j)
        j=j-1
      end do
      values(j+1)=value
    end do
  end subroutine sort_real_ascending

!###############################################################################

  subroutine solve_general_gmres(apply_operator,rhs,solution,tolerance, &
      max_iterations,restart,status)
    procedure(mrsf_tda_operator) :: apply_operator
    real(kind=dp), intent(in) :: rhs(:),tolerance
    real(kind=dp), intent(inout) :: solution(:)
    integer, intent(in) :: max_iterations,restart
    integer, intent(out) :: status

    real(kind=dp), allocatable :: basis(:,:),hessenberg(:,:),givens_c(:), &
      givens_s(:),g(:),y(:),work(:),applied(:),residual(:)
    real(kind=dp) :: beta,denominator,norm_rhs,temp
    integer :: i,j,k,m,n,iterations,operator_status,used
    logical :: converged

    n=size(rhs); m=min(restart,n); status=0; iterations=0; converged=.false.
    if(n<=0 .or. size(solution)/=n .or. m<=0 .or. &
       .not.ieee_is_finite(tolerance) .or. tolerance<=0.0_dp .or. &
       any(.not.ieee_is_finite(rhs)) .or. &
       any(.not.ieee_is_finite(solution))) then
      status=-1
      return
    end if
    allocate(basis(n,m+1),hessenberg(m+1,m),givens_c(m),givens_s(m), &
      g(m+1),y(m),work(n),applied(n),residual(n))
    norm_rhs=sqrt(dot_product(rhs,rhs))
    if(norm_rhs<=tolerance) then
      deallocate(basis,hessenberg,givens_c,givens_s,g,y,work,applied,residual)
      return
    end if
    do while(iterations<max_iterations .and. .not.converged)
      call apply_operator(solution,applied,operator_status)
      if(operator_status/=0 .or. any(.not.ieee_is_finite(applied))) then
        status=-2
        exit
      end if
      residual=rhs-applied
      beta=sqrt(dot_product(residual,residual))
      if(beta<=tolerance*max(1.0_dp,norm_rhs)) then
        converged=.true.
        exit
      end if
      basis=0.0_dp; hessenberg=0.0_dp; g=0.0_dp
      givens_c=0.0_dp; givens_s=0.0_dp
      basis(:,1)=residual/beta; g(1)=beta; used=0
      do j=1,m
        if(iterations>=max_iterations) exit
        iterations=iterations+1; used=j
        call apply_operator(basis(:,j),work,operator_status)
        if(operator_status/=0 .or. any(.not.ieee_is_finite(work))) then
          status=-2
          exit
        end if
        do k=1,2
          do i=1,j
            temp=dot_product(basis(:,i),work)
            hessenberg(i,j)=hessenberg(i,j)+temp
            work=work-temp*basis(:,i)
          end do
        end do
        hessenberg(j+1,j)=sqrt(dot_product(work,work))
        if(.not.ieee_is_finite(hessenberg(j+1,j))) then
          status=-6
          exit
        end if
        if(hessenberg(j+1,j)>100.0_dp*epsilon(1.0_dp)) &
          basis(:,j+1)=work/hessenberg(j+1,j)
        do i=1,j-1
          temp=givens_c(i)*hessenberg(i,j)+givens_s(i)*hessenberg(i+1,j)
          hessenberg(i+1,j)=-givens_s(i)*hessenberg(i,j)+ &
            givens_c(i)*hessenberg(i+1,j)
          hessenberg(i,j)=temp
        end do
        denominator=hypot(hessenberg(j,j),hessenberg(j+1,j))
        if(denominator<=100.0_dp*epsilon(1.0_dp)) then
          if(abs(g(j+1))<=tolerance*max(1.0_dp,norm_rhs)) then
            givens_c(j)=1.0_dp; givens_s(j)=0.0_dp
          else
            status=-3
          end if
        else
          givens_c(j)=hessenberg(j,j)/denominator
          givens_s(j)=hessenberg(j+1,j)/denominator
        end if
        if(status/=0) exit
        hessenberg(j,j)=givens_c(j)*hessenberg(j,j)+ &
          givens_s(j)*hessenberg(j+1,j)
        hessenberg(j+1,j)=0.0_dp
        temp=givens_c(j)*g(j)+givens_s(j)*g(j+1)
        g(j+1)=-givens_s(j)*g(j)+givens_c(j)*g(j+1)
        g(j)=temp
        if(abs(g(j+1))<=tolerance*max(1.0_dp,norm_rhs)) then
          converged=.true.
          exit
        end if
      end do
      if(status/=0 .or. used<=0) exit
      y(1:used)=g(1:used)
      do i=used,1,-1
        if(abs(hessenberg(i,i))<=100.0_dp*epsilon(1.0_dp)) then
          status=-4
          exit
        end if
        if(i<used) y(i)=y(i)-dot_product(hessenberg(i,i+1:used),y(i+1:used))
        y(i)=y(i)/hessenberg(i,i)
      end do
      if(status/=0) exit
      solution=solution+matmul(basis(:,1:used),y(1:used))
      if(any(.not.ieee_is_finite(solution))) then
        status=-6
        exit
      end if
    end do
    if(status==0) then
      call apply_operator(solution,applied,operator_status)
      if(operator_status/=0 .or. any(.not.ieee_is_finite(applied)) .or. &
         sqrt(dot_product(rhs-applied,rhs-applied))> &
           tolerance*max(1.0_dp,norm_rhs)) status=-5
    end if
    deallocate(basis,hessenberg,givens_c,givens_s,g,y,work,applied,residual)
  end subroutine solve_general_gmres

!###############################################################################

  subroutine solve_projected_gmres(apply_operator,omega,x0,rhs,solution, &
      tolerance,max_iterations,restart,status)
    procedure(mrsf_tda_operator) :: apply_operator
    real(kind=dp), intent(in) :: omega,x0(:),rhs(:),tolerance
    real(kind=dp), intent(inout) :: solution(:)
    integer, intent(in) :: max_iterations,restart
    integer, intent(out) :: status

    real(kind=dp), allocatable :: basis(:,:),hessenberg(:,:),givens_c(:), &
      givens_s(:),g(:),y(:),work(:),amat_work(:),residual(:)
    real(kind=dp) :: beta,denominator,norm_rhs,temp
    integer :: i,j,k,m,n,iterations,operator_status,used
    logical :: converged

    n=size(rhs); m=min(restart,n); status=0; iterations=0; converged=.false.
    if(n<=0 .or. size(x0)/=n .or. size(solution)/=n .or. m<=0 .or. &
       .not.ieee_is_finite(omega) .or. &
       .not.ieee_is_finite(tolerance) .or. tolerance<=0.0_dp .or. &
       any(.not.ieee_is_finite(x0)) .or. &
       any(.not.ieee_is_finite(rhs)) .or. &
       any(.not.ieee_is_finite(solution))) then
      status=-1
      return
    end if
    allocate(basis(n,m+1),hessenberg(m+1,m),givens_c(m),givens_s(m), &
      g(m+1),y(m),work(n),amat_work(n),residual(n))
    norm_rhs=sqrt(dot_product(rhs,rhs))
    if(.not.ieee_is_finite(norm_rhs)) then
      deallocate(basis,hessenberg,givens_c,givens_s,g,y,work,amat_work,residual)
      status=-6
      return
    end if
    if(norm_rhs<=tolerance) then
      deallocate(basis,hessenberg,givens_c,givens_s,g,y,work,amat_work,residual)
      return
    end if
    do while(iterations<max_iterations .and. .not.converged)
      call apply_projected_operator(apply_operator,omega,x0,solution, &
        amat_work,operator_status)
      if(operator_status/=0 .or. any(.not.ieee_is_finite(amat_work))) then
        status=-2
        exit
      end if
      residual=rhs-amat_work
      beta=sqrt(dot_product(residual,residual))
      if(.not.ieee_is_finite(beta)) then
        status=-6
        exit
      end if
      if(beta<=tolerance*max(1.0_dp,norm_rhs)) then
        converged=.true.
        exit
      end if
      basis=0.0_dp; hessenberg=0.0_dp; g=0.0_dp
      givens_c=0.0_dp; givens_s=0.0_dp
      basis(:,1)=residual/beta; g(1)=beta; used=0
      do j=1,m
        if(iterations>=max_iterations) exit
        iterations=iterations+1; used=j
        call apply_projected_operator(apply_operator,omega,x0,basis(:,j), &
          work,operator_status)
        if(operator_status/=0 .or. any(.not.ieee_is_finite(work))) then
          status=-2
          exit
        end if
        ! Two modified Gram--Schmidt passes keep the physical Krylov basis
        ! orthogonal when the target state is close to another root.
        do k=1,2
          do i=1,j
            temp=dot_product(basis(:,i),work)
            hessenberg(i,j)=hessenberg(i,j)+temp
            work=work-temp*basis(:,i)
          end do
        end do
        hessenberg(j+1,j)=sqrt(dot_product(work,work))
        if(.not.ieee_is_finite(hessenberg(j+1,j))) then
          status=-6
          exit
        end if
        if(hessenberg(j+1,j)>100.0_dp*epsilon(1.0_dp)) &
          basis(:,j+1)=work/hessenberg(j+1,j)
        do i=1,j-1
          temp=givens_c(i)*hessenberg(i,j)+givens_s(i)*hessenberg(i+1,j)
          hessenberg(i+1,j)=-givens_s(i)*hessenberg(i,j)+ &
            givens_c(i)*hessenberg(i+1,j)
          hessenberg(i,j)=temp
        end do
        denominator=hypot(hessenberg(j,j),hessenberg(j+1,j))
        if(.not.ieee_is_finite(denominator)) then
          status=-6
          exit
        end if
        if(denominator<=100.0_dp*epsilon(1.0_dp)) then
          if(abs(g(j+1))<=tolerance*max(1.0_dp,norm_rhs)) then
            givens_c(j)=1.0_dp; givens_s(j)=0.0_dp
          else
            status=-3
          end if
        else
          givens_c(j)=hessenberg(j,j)/denominator
          givens_s(j)=hessenberg(j+1,j)/denominator
        end if
        if(status/=0) exit
        hessenberg(j,j)=givens_c(j)*hessenberg(j,j)+ &
          givens_s(j)*hessenberg(j+1,j)
        hessenberg(j+1,j)=0.0_dp
        temp=givens_c(j)*g(j)+givens_s(j)*g(j+1)
        g(j+1)=-givens_s(j)*g(j)+givens_c(j)*g(j+1)
        g(j)=temp
        if(abs(g(j+1))<=tolerance*max(1.0_dp,norm_rhs)) then
          converged=.true.
          exit
        end if
      end do
      if(status/=0 .or. used<=0) exit
      y(1:used)=g(1:used)
      do i=used,1,-1
        if(abs(hessenberg(i,i))<=100.0_dp*epsilon(1.0_dp)) then
          status=-4
          exit
        end if
        if(i<used) y(i)=y(i)-dot_product(hessenberg(i,i+1:used),y(i+1:used))
        y(i)=y(i)/hessenberg(i,i)
        if(.not.ieee_is_finite(y(i))) then
          status=-6
          exit
        end if
      end do
      if(status/=0) exit
      solution=solution+matmul(basis(:,1:used),y(1:used))
      solution=solution-x0*dot_product(x0,solution)
      if(any(.not.ieee_is_finite(solution))) then
        status=-6
        exit
      end if
    end do
    if(status==0) then
      call apply_projected_operator(apply_operator,omega,x0,solution, &
        amat_work,operator_status)
      if(operator_status/=0 .or. any(.not.ieee_is_finite(amat_work)) .or. &
         sqrt(dot_product(rhs-amat_work,rhs-amat_work))> &
           tolerance*max(1.0_dp,norm_rhs)) status=-5
    end if
    deallocate(basis,hessenberg,givens_c,givens_s,g,y,work,amat_work,residual)
  end subroutine solve_projected_gmres

!###############################################################################

  subroutine apply_projected_operator(apply_operator,omega,x0,vector,result,status)
    procedure(mrsf_tda_operator) :: apply_operator
    real(kind=dp), intent(in) :: omega,x0(:),vector(:)
    real(kind=dp), intent(out) :: result(:)
    integer, intent(out) :: status
    real(kind=dp), allocatable :: projected(:),applied(:)
    real(kind=dp) :: parallel_component

    allocate(projected(size(vector)),applied(size(vector)))
    status=0
    result=0.0_dp
    if(size(x0)/=size(vector) .or. size(result)/=size(vector) .or. &
       .not.ieee_is_finite(omega) .or. any(.not.ieee_is_finite(x0)) .or. &
       any(.not.ieee_is_finite(vector))) then
      status=-1
      deallocate(projected,applied)
      return
    end if
    parallel_component=dot_product(x0,vector)
    projected=vector-parallel_component*x0
    call apply_operator(projected,applied,status)
    if(status==0 .and. all(ieee_is_finite(applied))) then
      result=applied-omega*projected
      result=result-x0*dot_product(x0,result)+parallel_component*x0
      if(any(.not.ieee_is_finite(result))) status=-1
    else
      status=-1
      result=0.0_dp
    end if
    deallocate(projected,applied)
  end subroutine apply_projected_operator

!###############################################################################

  subroutine solve_linear_pivot(matrix, rhs, solution, threshold, status)
    ! Partial-pivoting Gaussian elimination used only by the dense reference
    ! response above.  Keeping it local avoids linking a production BLAS/LAPACK
    ! implementation into the standalone algebraic verification program.

    real(kind=dp), intent(in) :: matrix(:,:), rhs(:), threshold
    real(kind=dp), intent(out) :: solution(:)
    integer, intent(out) :: status

    real(kind=dp), allocatable :: work(:,:), vector(:), row(:)
    real(kind=dp) :: factor, scale
    integer :: i, k, n, pivot

    n = size(rhs)
    status = 0
    solution = 0.0_dp
    if (n <= 0 .or. any(shape(matrix) /= [n,n]) .or. &
        size(solution) /= n) then
      status = -1
      return
    end if
    allocate(work(n,n),vector(n),row(n))
    work = matrix
    vector = rhs
    scale = max(1.0_dp,maxval(abs(work)))
    do k = 1,n-1
      pivot = k-1+maxloc(abs(work(k:n,k)),dim=1)
      if (abs(work(pivot,k)) <= threshold*scale) then
        status = k
        exit
      end if
      if (pivot /= k) then
        row = work(k,:)
        work(k,:) = work(pivot,:)
        work(pivot,:) = row
        factor = vector(k)
        vector(k) = vector(pivot)
        vector(pivot) = factor
      end if
      do i = k+1,n
        factor = work(i,k)/work(k,k)
        work(i,k:n) = work(i,k:n)-factor*work(k,k:n)
        vector(i) = vector(i)-factor*vector(k)
      end do
    end do
    if (status == 0 .and. abs(work(n,n)) <= threshold*scale) status = n
    if (status == 0) then
      do i = n,1,-1
        solution(i) = (vector(i)-dot_product(work(i,i+1:n), &
          solution(i+1:n)))/work(i,i)
      end do
    end if
    deallocate(work,vector,row)
  end subroutine solve_linear_pivot

!###############################################################################

  subroutine complete_rhf_orbital_response(mo, nocc, sx_mo, uov, umat, dmo_occ, dp_ao, status)
    ! Complete the occupied-virtual CPHF amplitudes into a full first-order
    ! MO connection.  The occupied-occupied and virtual-virtual rotations use
    ! the symmetric gauge, U_pq=-S^x_pq/2.  The transposed occupied-virtual
    ! block follows from differentiated orthonormality,
    !
    !                 U + U^T + S^x = 0.
    !
    ! The returned AO density is the spin-summed closed-shell response.

    real(kind=dp), intent(in) :: mo(:,:), sx_mo(:,:), uov(:)
    integer, intent(in) :: nocc
    real(kind=dp), intent(out) :: umat(:,:), dmo_occ(:,:), dp_ao(:,:)
    integer, intent(out) :: status

    integer :: a, i, ia, nbf, nvir

    nbf = size(mo, 1)
    nvir = nbf - nocc
    status = 0
    umat = 0.0_dp
    dmo_occ = 0.0_dp
    dp_ao = 0.0_dp

    if (nbf <= 0 .or. size(mo,2) /= nbf .or. nocc <= 0 .or. nocc >= nbf .or. &
        any(shape(sx_mo) /= [nbf,nbf]) .or. size(uov) /= nocc*nvir .or. &
        any(shape(umat) /= [nbf,nbf]) .or. &
        any(shape(dmo_occ) /= [nbf,nocc]) .or. &
        any(shape(dp_ao) /= [nbf,nbf])) then
      status = -1
      return
    end if

    umat(1:nocc,1:nocc) = -0.5_dp*sx_mo(1:nocc,1:nocc)
    umat(nocc+1:nbf,nocc+1:nbf) = &
      -0.5_dp*sx_mo(nocc+1:nbf,nocc+1:nbf)
    ia = 0
    do a = 1, nvir
      do i = 1, nocc
        ia = ia + 1
        umat(nocc+a,i) = uov(ia)
        umat(i,nocc+a) = -sx_mo(i,nocc+a) - uov(ia)
      end do
    end do

    dmo_occ = matmul(mo, umat(:,1:nocc))
    dp_ao = 2.0_dp*(matmul(dmo_occ, transpose(mo(:,1:nocc))) + &
                    matmul(mo(:,1:nocc), transpose(dmo_occ)))
  end subroutine complete_rhf_orbital_response

!###############################################################################

  subroutine assemble_tdhf_sigma_derivative(z, nocc, umat, eps_deriv, gmat, &
                                             deri_ov, inner_ov, sigma_deriv, status)
    ! Assemble one differentiated TDHF/TDDFT response-operator action.
    ! For either the A-B or A+B channel and every nuclear coordinate K,
    !
    ! dSigma_ia = dERI_ia + sum_p(U_pi G_pa + G_ip U_pa)
    !             + W[dP]_ia + (eps_a^K-eps_i^K) Z_ia.
    !
    ! The caller supplies channel-specific derivative-ERI, G, and inner
    ! density-response contractions.  Keeping this contraction independent
    ! of their integral realization makes every primitive directly testable.

    real(kind=dp), intent(in) :: z(:), umat(:,:,:), eps_deriv(:,:), gmat(:,:)
    integer, intent(in) :: nocc
    real(kind=dp), intent(in) :: deri_ov(:,:), inner_ov(:,:)
    real(kind=dp), intent(out) :: sigma_deriv(:,:)
    integer, intent(out) :: status

    integer :: a, i, ia, k, ncoord, nexc, nmo, nvir

    nmo = size(gmat,1)
    ncoord = size(umat,3)
    nexc = size(z)
    status = 0
    sigma_deriv = 0.0_dp
    if (nmo <= 1 .or. nocc <= 0 .or. nocc >= nmo .or. &
        size(gmat,2) /= nmo .or. size(umat,1) /= nmo .or. &
        size(umat,2) /= nmo .or. any(shape(eps_deriv) /= [nmo,ncoord]) .or. &
        any(shape(deri_ov) /= [nexc,ncoord]) .or. &
        any(shape(inner_ov) /= [nexc,ncoord]) .or. &
        any(shape(sigma_deriv) /= [nexc,ncoord])) then
      status = -1
      return
    end if

    nvir = nmo - nocc
    if (nexc /= nocc*nvir) then
      status = -2
      return
    end if

    ia = 0
    do a = 1, nvir
      do i = 1, nocc
        ia = ia + 1
        do k = 1, ncoord
          sigma_deriv(ia,k) = deri_ov(ia,k) + inner_ov(ia,k) &
            + dot_product(umat(:,i,k), gmat(:,nocc+a)) &
            + dot_product(gmat(i,:), umat(:,nocc+a,k)) &
            + (eps_deriv(nocc+a,k)-eps_deriv(i,k))*z(ia)
        end do
      end do
    end do
  end subroutine assemble_tdhf_sigma_derivative

!###############################################################################

  subroutine solve_tdhf_amplitude_response(amb, apb, omega, u0, v0, &
                                            dambu, dapbv, du, dv, domega, &
                                            residual_max, status, tol, maxit)
    ! Solve the differentiated full-response TDHF/TDDFT equations
    !
    !   (A-B) dU - omega dV = domega V - d(A-B) U,
    !   (A+B) dV - omega dU = domega U - d(A+B) V.
    !
    ! The coupled matrix has the null vector (U,V).  Each right-hand side
    ! is projected perpendicular to that vector.  After the projected
    ! conjugate-gradient solution, the remaining null-vector component is
    ! chosen so that d(U^T V)=0.  This is the same gauge condition used by
    ! the GAMESS TDHBLD implementation.

    real(kind=dp), intent(in) :: amb(:,:), apb(:,:)
    real(kind=dp), intent(in) :: omega, u0(:), v0(:)
    real(kind=dp), intent(in) :: dambu(:,:), dapbv(:,:)
    real(kind=dp), intent(out) :: du(:,:), dv(:,:), domega(:)
    real(kind=dp), intent(out) :: residual_max
    integer, intent(out) :: status
    real(kind=dp), intent(in), optional :: tol
    integer, intent(in), optional :: maxit

    real(kind=dp), allocatable :: b(:), x(:), r(:), p(:), q(:), nullvec(:)
    real(kind=dp) :: alpha, beta, bnorm, gauge_denominator, gauge_value
    real(kind=dp) :: n2, pq, rr, rr_new, solve_tol, uvnorm
    integer :: i, iteration, ncoord, nexc, niter

    nexc = size(u0)
    ncoord = size(dambu, 2)
    status = 0
    residual_max = 0.0_dp
    domega = 0.0_dp
    du = 0.0_dp
    dv = 0.0_dp

    if (nexc <= 0 .or. size(v0) /= nexc .or. &
        any(shape(amb) /= [nexc, nexc]) .or. &
        any(shape(apb) /= [nexc, nexc]) .or. &
        size(dambu, 1) /= nexc .or. &
        any(shape(dapbv) /= [nexc, ncoord]) .or. &
        any(shape(du) /= [nexc, ncoord]) .or. &
        any(shape(dv) /= [nexc, ncoord]) .or. &
        size(domega) /= ncoord) then
      status = -1
      return
    end if

    solve_tol = 1.0e-12_dp
    if (present(tol)) solve_tol = tol
    niter = 400
    if (present(maxit)) niter = maxit
    if (solve_tol <= 0.0_dp .or. niter <= 0) then
      status = -2
      return
    end if

    uvnorm = dot_product(u0, v0)
    gauge_denominator = 2.0_dp*uvnorm
    if (abs(gauge_denominator) <= tiny(1.0_dp)) then
      status = -3
      return
    end if

    allocate(b(2*nexc), x(2*nexc), r(2*nexc), p(2*nexc), q(2*nexc), &
             nullvec(2*nexc))
    nullvec(1:nexc) = u0
    nullvec(nexc+1:2*nexc) = v0
    n2 = dot_product(nullvec, nullvec)

    do i = 1, ncoord
      domega(i) = (dot_product(u0, dambu(:,i)) + &
                   dot_product(v0, dapbv(:,i)))/(2.0_dp*uvnorm)
      b(1:nexc) = domega(i)*v0 - dambu(:,i)
      b(nexc+1:2*nexc) = domega(i)*u0 - dapbv(:,i)
      call project_from_null(b, nullvec, n2)

      bnorm = sqrt(dot_product(b, b))
      if (bnorm <= tiny(1.0_dp)) cycle

      x = 0.0_dp
      r = b
      p = r
      rr = dot_product(r, r)

      do iteration = 1, niter
        call apply_coupled_matrix(amb, apb, omega, p, q, nexc)
        call project_from_null(q, nullvec, n2)
        pq = dot_product(p, q)
        if (pq <= tiny(1.0_dp)) then
          status = i
          exit
        end if
        alpha = rr/pq
        x = x + alpha*p
        r = r - alpha*q
        rr_new = dot_product(r, r)
        if (sqrt(rr_new) <= solve_tol*bnorm) then
          rr = rr_new
          exit
        end if
        beta = rr_new/rr
        p = r + beta*p
        rr = rr_new
      end do
      if (status /= 0) exit
      if (sqrt(rr) > solve_tol*bnorm) then
        status = i
        exit
      end if

      ! Enforce V^T dU + U^T dV = 0 without changing the equations.
      gauge_value = dot_product(v0, x(1:nexc)) + &
                    dot_product(u0, x(nexc+1:2*nexc))
      x = x - (gauge_value/gauge_denominator)*nullvec
      du(:,i) = x(1:nexc)
      dv(:,i) = x(nexc+1:2*nexc)

      call apply_coupled_matrix(amb, apb, omega, x, q, nexc)
      residual_max = max(residual_max, maxval(abs(q-b)))
    end do

    deallocate(b, x, r, p, q, nullvec)

  end subroutine solve_tdhf_amplitude_response

!###############################################################################

  subroutine solve_tdhf_z_response(orbital_hessian, rhs_derivative, &
                                    operator_derivative_z, dz, residual_max, &
                                    status, tol, maxit)
    ! Differentiate H Z = R for every Cartesian perturbation:
    !
    !   H dZ = dR - (dH) Z.

    real(kind=dp), intent(in) :: orbital_hessian(:,:)
    real(kind=dp), intent(in) :: rhs_derivative(:,:), operator_derivative_z(:,:)
    real(kind=dp), intent(out) :: dz(:,:)
    real(kind=dp), intent(out) :: residual_max
    integer, intent(out) :: status
    real(kind=dp), intent(in), optional :: tol
    integer, intent(in), optional :: maxit

    real(kind=dp), allocatable :: b(:), p(:), q(:), r(:), x(:)
    real(kind=dp) :: alpha, beta, bnorm, pq, rr, rr_new, solve_tol
    integer :: i, iteration, ncoord, nexc, niter

    nexc = size(orbital_hessian, 1)
    ncoord = size(rhs_derivative, 2)
    status = 0
    residual_max = 0.0_dp
    dz = 0.0_dp

    if (nexc <= 0 .or. size(orbital_hessian, 2) /= nexc .or. &
        size(rhs_derivative, 1) /= nexc .or. &
        any(shape(operator_derivative_z) /= [nexc, ncoord]) .or. &
        any(shape(dz) /= [nexc, ncoord])) then
      status = -1
      return
    end if

    solve_tol = 1.0e-12_dp
    if (present(tol)) solve_tol = tol
    niter = 400
    if (present(maxit)) niter = maxit
    if (solve_tol <= 0.0_dp .or. niter <= 0) then
      status = -2
      return
    end if

    allocate(b(nexc), p(nexc), q(nexc), r(nexc), x(nexc))
    do i = 1, ncoord
      b = rhs_derivative(:,i) - operator_derivative_z(:,i)
      bnorm = sqrt(dot_product(b, b))
      if (bnorm <= tiny(1.0_dp)) cycle
      x = 0.0_dp
      r = b
      p = r
      rr = dot_product(r, r)

      do iteration = 1, niter
        q = matmul(orbital_hessian, p)
        pq = dot_product(p, q)
        if (pq <= tiny(1.0_dp)) then
          status = i
          exit
        end if
        alpha = rr/pq
        x = x + alpha*p
        r = r - alpha*q
        rr_new = dot_product(r, r)
        if (sqrt(rr_new) <= solve_tol*bnorm) then
          rr = rr_new
          exit
        end if
        beta = rr_new/rr
        p = r + beta*p
        rr = rr_new
      end do
      if (status /= 0) exit
      if (sqrt(rr) > solve_tol*bnorm) then
        status = i
        exit
      end if
      dz(:,i) = x
      residual_max = max(residual_max, &
        maxval(abs(matmul(orbital_hessian, x)-b)))
    end do
    deallocate(b, p, q, r, x)

  end subroutine solve_tdhf_z_response

!###############################################################################

  subroutine solve_mrsf_z_response_matrix_free(apply_orbital_hessian, &
      rhs_derivative,operator_derivative_z,dz,residual_max,status,tol,maxit, &
      restart)
    ! Differentiate the spin-adapted MRSF orbital Lagrange-multiplier equation,
    !
    !   H_orb dZ(K) = dR(K) - dH_orb(K) Z.
    !
    ! The ROHF/ROKS orbital Hessian need not be positive definite.  Restarted
    ! GMRES therefore replaces the conjugate-gradient assumption of the older
    ! dense TDHF reference routine.  The production caller supplies the same
    ! matrix-free orbital-Hessian action used by the gradient Z-vector.

    procedure(mrsf_tda_operator) :: apply_orbital_hessian
    real(kind=dp), intent(in) :: rhs_derivative(:,:),operator_derivative_z(:,:)
    real(kind=dp), intent(out) :: dz(:,:),residual_max
    integer, intent(out) :: status
    real(kind=dp), intent(in), optional :: tol
    integer, intent(in), optional :: maxit,restart

    real(kind=dp), allocatable :: rhs(:),solution(:),applied(:)
    real(kind=dp) :: solve_tol
    integer :: coordinate,n,ncoord,niter,nrestart,operator_status,solve_status

    n=size(rhs_derivative,1)
    ncoord=size(rhs_derivative,2)
    solve_tol=1.0e-11_dp
    if(present(tol)) solve_tol=tol
    niter=max(200,4*n)
    if(present(maxit)) niter=maxit
    nrestart=min(max(12,n/4),max(1,n))
    if(present(restart)) nrestart=restart
    status=0
    residual_max=0.0_dp
    dz=0.0_dp
    if(n<=0 .or. ncoord<=0 .or. &
       any(shape(operator_derivative_z)/=[n,ncoord]) .or. &
       any(shape(dz)/=[n,ncoord]) .or. solve_tol<=0.0_dp .or. &
       niter<=0 .or. nrestart<=0 .or. &
       any(.not.ieee_is_finite(rhs_derivative)) .or. &
       any(.not.ieee_is_finite(operator_derivative_z))) then
      status=-1
      return
    end if

    allocate(rhs(n),solution(n),applied(n))
    do coordinate=1,ncoord
      rhs=rhs_derivative(:,coordinate)-operator_derivative_z(:,coordinate)
      solution=0.0_dp
      call solve_general_gmres(apply_orbital_hessian,rhs,solution,solve_tol, &
        niter,nrestart,solve_status)
      if(solve_status/=0) then
        write(error_unit,'(A,I0,A,I0,A,1P,E12.4)') &
          'MRSF differentiated Z solve failed for coordinate ',coordinate, &
          ' with status ',solve_status,', max|rhs|=',maxval(abs(rhs))
        status=coordinate
        exit
      end if
      dz(:,coordinate)=solution
      call apply_orbital_hessian(solution,applied,operator_status)
      if(operator_status/=0 .or. any(.not.ieee_is_finite(applied)) .or. &
         any(.not.ieee_is_finite(solution))) then
        status=coordinate
        exit
      end if
      residual_max=max(residual_max,maxval(abs(applied-rhs)))
    end do
    deallocate(rhs,solution,applied)
  end subroutine solve_mrsf_z_response_matrix_free

!###############################################################################

  subroutine solve_mrsf_z_response_batch_matrix_free(apply_orbital_hessian, &
      rhs_derivative,operator_derivative_z,dz,residual_max,status,tol,maxit, &
      apply_preconditioner)
    ! Batched differentiated MRSF orbital Z-vector equations,
    !
    !   H_orb dZ(K) = dR(K) - dH_orb(K) Z.
    !
    ! All active nuclear perturbations share one AO integral/XC response pass
    ! per simultaneous Krylov step.  This is the multi-RHS analogue of the
    ! GAMESS AO-CPKS organization and never forms the orbital Hessian.

    procedure(mrsf_tda_batch_operator) :: apply_orbital_hessian
    procedure(mrsf_tda_batch_operator), optional :: apply_preconditioner
    real(kind=dp), intent(in) :: rhs_derivative(:,:),operator_derivative_z(:,:)
    real(kind=dp), intent(out) :: dz(:,:),residual_max
    integer, intent(out) :: status
    real(kind=dp), intent(in), optional :: tol
    integer, intent(in), optional :: maxit

    real(kind=dp), allocatable :: rhs(:,:)
    real(kind=dp) :: solve_tol
    integer :: n,ncoord,niter

    n=size(rhs_derivative,1)
    ncoord=size(rhs_derivative,2)
    solve_tol=1.0e-10_dp
    if(present(tol)) solve_tol=tol
    niter=max(200,4*n)
    if(present(maxit)) niter=maxit
    status=0
    residual_max=0.0_dp
    dz=0.0_dp
    if(n<=0 .or. ncoord<=0 .or. &
       any(shape(operator_derivative_z)/=[n,ncoord]) .or. &
       any(shape(dz)/=[n,ncoord]) .or. &
       .not.ieee_is_finite(solve_tol) .or. solve_tol<=0.0_dp .or. &
       niter<=0 .or. any(.not.ieee_is_finite(rhs_derivative)) .or. &
       any(.not.ieee_is_finite(operator_derivative_z))) then
      status=-1
      return
    end if
    allocate(rhs(n,ncoord))
    rhs=rhs_derivative-operator_derivative_z
    if(present(apply_preconditioner)) then
      ! The MRSF ROHF orbital Hessian is generally indefinite, so a PCG trial
      ! can break down after spending several exact AO builds and must then
      ! restart from zero.  Begin directly with simultaneous GMRES.  The
      ! diagonal preconditioner changes convergence only; the returned vectors
      ! are certified against the unpreconditioned equations.
      call solve_batched_gmres(apply_orbital_hessian,rhs,dz, &
        residual_max,status,solve_tol,niter, &
        apply_preconditioner=apply_preconditioner)
    else
      call solve_batched_gmres(apply_orbital_hessian,rhs,dz, &
        residual_max,status,solve_tol,niter)
    end if
    deallocate(rhs)
  end subroutine solve_mrsf_z_response_batch_matrix_free

!###############################################################################

  subroutine solve_batched_pcg(apply_operator,apply_preconditioner,rhs, &
      solution,residual_max,status,tolerance,max_iterations)
    ! Low-storage simultaneous PCG for a symmetric orbital Hessian.  Each
    ! perturbation retains its own scalar recurrence, but all unconverged
    ! search directions are packed into one exact AO operator application.
    ! This is the multiple-response organization of GAMESS AOCPCG, generalized
    ! to the MRSF ROHF Z-vector space and followed by an exact batch residual
    ! certification.  The preconditioner affects convergence only.
    procedure(mrsf_tda_batch_operator) :: apply_operator
    procedure(mrsf_tda_batch_operator) :: apply_preconditioner
    real(kind=dp), intent(in) :: rhs(:,:),tolerance
    real(kind=dp), intent(out) :: solution(:,:),residual_max
    integer, intent(out) :: status
    integer, intent(in) :: max_iterations

    real(kind=dp), allocatable :: residual(:,:),preconditioned(:,:), &
      search(:,:),applied(:,:),active_vectors(:,:),active_images(:,:), &
      norm_rhs(:),rho(:),rho_new(:)
    integer, allocatable :: active_columns(:)
    logical, allocatable :: converged(:)
    real(kind=dp) :: alpha,beta,component_residual,denominator,current_norm,scale
    integer :: active_index,coordinate,iteration,n,nactive,ncoord,operator_status

    n=size(rhs,1)
    ncoord=size(rhs,2)
    status=0
    residual_max=0.0_dp
    solution=0.0_dp
    if(n<=0 .or. ncoord<=0 .or. any(shape(solution)/=[n,ncoord]) .or. &
       .not.ieee_is_finite(tolerance) .or. tolerance<=0.0_dp .or. &
       max_iterations<=0 .or. any(.not.ieee_is_finite(rhs))) then
      status=-1
      return
    end if
    allocate(residual(n,ncoord),preconditioned(n,ncoord), &
      search(n,ncoord),applied(n,ncoord),active_vectors(n,ncoord), &
      active_images(n,ncoord),norm_rhs(ncoord),rho(ncoord), &
      rho_new(ncoord),active_columns(ncoord),converged(ncoord))
    do coordinate=1,ncoord
      norm_rhs(coordinate)=sqrt(dot_product(rhs(:,coordinate), &
        rhs(:,coordinate)))
    end do
    if(any(.not.ieee_is_finite(norm_rhs))) then
      status=-2
      go to 900
    end if
    converged=norm_rhs<=tolerance
    if(all(converged)) go to 800

    ! Diagonal-response initial guess, as in ordinary AO CPHF.
    call apply_preconditioner(rhs,solution,operator_status)
    if(operator_status/=0 .or. any(.not.ieee_is_finite(solution))) then
      status=-3
      go to 900
    end if
    call apply_operator(solution,applied,operator_status)
    if(operator_status/=0 .or. any(.not.ieee_is_finite(applied))) then
      status=-4
      go to 900
    end if
    residual=rhs-applied
    do coordinate=1,ncoord
      current_norm=sqrt(dot_product(residual(:,coordinate), &
        residual(:,coordinate)))
      converged(coordinate)=current_norm<= &
        tolerance*max(1.0_dp,norm_rhs(coordinate))
    end do
    if(all(converged)) go to 800
    call apply_preconditioner(residual,preconditioned,operator_status)
    if(operator_status/=0 .or. any(.not.ieee_is_finite(preconditioned))) then
      status=-3
      go to 900
    end if
    search=preconditioned
    do coordinate=1,ncoord
      rho(coordinate)=dot_product(residual(:,coordinate), &
        preconditioned(:,coordinate))
    end do

    do iteration=1,max_iterations
      nactive=0
      do coordinate=1,ncoord
        if(.not.converged(coordinate)) then
          nactive=nactive+1
          active_columns(nactive)=coordinate
          active_vectors(:,nactive)=search(:,coordinate)
        end if
      end do
      if(nactive==0) exit
      call apply_operator(active_vectors(:,1:nactive), &
        active_images(:,1:nactive),operator_status)
      if(operator_status/=0 .or. &
         any(.not.ieee_is_finite(active_images(:,1:nactive)))) then
        status=-4
        go to 900
      end if
      do active_index=1,nactive
        coordinate=active_columns(active_index)
        applied(:,coordinate)=active_images(:,active_index)
        denominator=dot_product(search(:,coordinate),applied(:,coordinate))
        scale=max(1.0_dp,sqrt(dot_product(search(:,coordinate), &
          search(:,coordinate)))*sqrt(dot_product(applied(:,coordinate), &
          applied(:,coordinate))))
        if(.not.ieee_is_finite(rho(coordinate)) .or. &
           .not.ieee_is_finite(denominator) .or. &
           abs(rho(coordinate))<=1.0e-35_dp .or. &
           abs(denominator)<=100.0_dp*epsilon(1.0_dp)*scale) then
          status=coordinate
          go to 900
        end if
        alpha=rho(coordinate)/denominator
        solution(:,coordinate)=solution(:,coordinate)+ &
          alpha*search(:,coordinate)
        residual(:,coordinate)=residual(:,coordinate)- &
          alpha*applied(:,coordinate)
        current_norm=sqrt(dot_product(residual(:,coordinate), &
          residual(:,coordinate)))
        if(.not.ieee_is_finite(current_norm)) then
          status=coordinate
          go to 900
        end if
        converged(coordinate)=current_norm<= &
          tolerance*max(1.0_dp,norm_rhs(coordinate))
      end do
      if(all(converged)) exit

      nactive=0
      do coordinate=1,ncoord
        if(.not.converged(coordinate)) then
          nactive=nactive+1
          active_columns(nactive)=coordinate
          active_vectors(:,nactive)=residual(:,coordinate)
        end if
      end do
      call apply_preconditioner(active_vectors(:,1:nactive), &
        active_images(:,1:nactive),operator_status)
      if(operator_status/=0 .or. &
         any(.not.ieee_is_finite(active_images(:,1:nactive)))) then
        status=-3
        go to 900
      end if
      do active_index=1,nactive
        coordinate=active_columns(active_index)
        preconditioned(:,coordinate)=active_images(:,active_index)
        rho_new(coordinate)=dot_product(residual(:,coordinate), &
          preconditioned(:,coordinate))
        if(.not.ieee_is_finite(rho_new(coordinate)) .or. &
           abs(rho(coordinate))<=1.0e-35_dp) then
          status=coordinate
          go to 900
        end if
        beta=rho_new(coordinate)/rho(coordinate)
        search(:,coordinate)=preconditioned(:,coordinate)+ &
          beta*search(:,coordinate)
        rho(coordinate)=rho_new(coordinate)
      end do
    end do
    if(.not.all(converged)) then
      do coordinate=1,ncoord
        if(.not.converged(coordinate)) then
          status=coordinate
          exit
        end if
      end do
      go to 900
    end if

800 continue
    call apply_operator(solution,applied,operator_status)
    if(operator_status/=0 .or. any(.not.ieee_is_finite(applied))) then
      status=-4
      go to 900
    end if
    do coordinate=1,ncoord
      current_norm=sqrt(dot_product(applied(:,coordinate)-rhs(:,coordinate), &
        applied(:,coordinate)-rhs(:,coordinate)))
      component_residual=maxval(abs(applied(:,coordinate)-rhs(:,coordinate)))
      residual_max=max(residual_max,component_residual)
      if(.not.ieee_is_finite(current_norm) .or. &
         component_residual>max(1.0e-8_dp,100.0_dp*tolerance)* &
           max(1.0_dp,norm_rhs(coordinate))) then
        status=coordinate
        go to 900
      end if
    end do

900 continue
    deallocate(residual,preconditioned,search,applied,active_vectors, &
      active_images,norm_rhs,rho,rho_new,active_columns,converged)
  end subroutine solve_batched_pcg

!###############################################################################

  subroutine solve_batched_gmres(apply_operator,rhs,solution,residual_max, &
      status,tolerance,max_iterations,apply_preconditioner)
    ! Simultaneous restarted right-preconditioned GMRES.  Each right-hand side
    ! keeps its own orthonormal Arnoldi recurrence, avoiding the conditioning
    ! loss that can occur when unrelated nuclear perturbations are forced into
    ! one least-squares problem.  At each Arnoldi level, however, all active
    ! directions are packed into one exact AO operator application.  Thus this
    ! retains the GAMESS multi-response organization without forming H or
    ! solving one nuclear coordinate at a time.
    procedure(mrsf_tda_batch_operator) :: apply_operator
    procedure(mrsf_tda_batch_operator), optional :: apply_preconditioner
    real(kind=dp), intent(in) :: rhs(:,:),tolerance
    real(kind=dp), intent(out) :: solution(:,:),residual_max
    integer, intent(out) :: status
    integer, intent(in) :: max_iterations

    real(kind=dp), allocatable :: basis(:,:,:),preconditioned_basis(:,:,:), &
      hessenberg(:,:,:),givens_c(:,:),givens_s(:,:),g(:,:),y(:,:), &
      residual(:,:),applied(:,:),active_vectors(:,:),active_preconditioned(:,:), &
      active_images(:,:),work(:),norm_rhs(:)
    integer, allocatable :: active_columns(:),used(:)
    logical, allocatable :: converged(:),cycle_done(:)
    real(kind=dp) :: beta,component_residual,current_norm,denominator, &
      orthogonal_norm,residual_scale,temp
    integer :: active_index,coordinate,cycle_start,i,iteration,j,k,m,n, &
      nactive,ncoord,operator_status,pass

    n=size(rhs,1)
    ncoord=size(rhs,2)
    status=0
    residual_max=0.0_dp
    solution=0.0_dp
    if(n<=0 .or. ncoord<=0 .or. any(shape(solution)/=[n,ncoord]) .or. &
       .not.ieee_is_finite(tolerance) .or. tolerance<=0.0_dp .or. &
       max_iterations<=0 .or. any(.not.ieee_is_finite(rhs))) then
      status=-1
      return
    end if

    m=min(n,min(64,max_iterations))
    allocate(basis(n,m+1,ncoord),preconditioned_basis(n,m,ncoord), &
      hessenberg(m+1,m,ncoord),givens_c(m,ncoord), &
      givens_s(m,ncoord),g(m+1,ncoord),y(m,ncoord), &
      residual(n,ncoord),applied(n,ncoord),active_vectors(n,ncoord), &
      active_preconditioned(n,ncoord),active_images(n,ncoord),work(n), &
      norm_rhs(ncoord),active_columns(ncoord),used(ncoord), &
      converged(ncoord),cycle_done(ncoord))
    do coordinate=1,ncoord
      norm_rhs(coordinate)=sqrt(dot_product(rhs(:,coordinate), &
        rhs(:,coordinate)))
    end do
    if(any(.not.ieee_is_finite(norm_rhs))) then
      status=-2
      go to 900
    end if

    iteration=0
    converged=norm_rhs<=tolerance
    do while(iteration<max_iterations .and. .not.all(converged))
      call apply_operator(solution,applied,operator_status)
      if(operator_status/=0 .or. any(.not.ieee_is_finite(applied))) then
        status=-3
        go to 900
      end if
      residual=rhs-applied
      do coordinate=1,ncoord
        current_norm=sqrt(dot_product(residual(:,coordinate), &
          residual(:,coordinate)))
        converged(coordinate)=current_norm<= &
          tolerance*max(1.0_dp,norm_rhs(coordinate))
      end do
      if(all(converged)) exit

      basis=0.0_dp
      preconditioned_basis=0.0_dp
      hessenberg=0.0_dp
      givens_c=0.0_dp
      givens_s=0.0_dp
      g=0.0_dp
      y=0.0_dp
      used=0
      cycle_done=converged
      do coordinate=1,ncoord
        if(converged(coordinate)) cycle
        beta=sqrt(dot_product(residual(:,coordinate),residual(:,coordinate)))
        if(.not.ieee_is_finite(beta) .or. beta<=tiny(1.0_dp)) then
          status=coordinate
          go to 900
        end if
        basis(:,1,coordinate)=residual(:,coordinate)/beta
        g(1,coordinate)=beta
      end do
      cycle_start=iteration

      do j=1,m
        if(iteration>=max_iterations .or. all(cycle_done)) exit
        nactive=0
        do coordinate=1,ncoord
          if(.not.cycle_done(coordinate)) then
            nactive=nactive+1
            active_columns(nactive)=coordinate
            active_vectors(:,nactive)=basis(:,j,coordinate)
          end if
        end do
        if(nactive<=0) exit
        if(present(apply_preconditioner)) then
          call apply_preconditioner(active_vectors(:,1:nactive), &
            active_preconditioned(:,1:nactive),operator_status)
          if(operator_status/=0 .or. &
             any(.not.ieee_is_finite(active_preconditioned(:,1:nactive)))) then
            status=-4
            go to 900
          end if
        else
          active_preconditioned(:,1:nactive)=active_vectors(:,1:nactive)
        end if
        call apply_operator(active_preconditioned(:,1:nactive), &
          active_images(:,1:nactive),operator_status)
        if(operator_status/=0 .or. &
           any(.not.ieee_is_finite(active_images(:,1:nactive)))) then
          status=-3
          go to 900
        end if

        do active_index=1,nactive
          coordinate=active_columns(active_index)
          preconditioned_basis(:,j,coordinate)= &
            active_preconditioned(:,active_index)
          work=active_images(:,active_index)
          do pass=1,2
            do k=1,j
              temp=dot_product(basis(:,k,coordinate),work)
              hessenberg(k,j,coordinate)=hessenberg(k,j,coordinate)+temp
              work=work-temp*basis(:,k,coordinate)
            end do
          end do
          orthogonal_norm=sqrt(dot_product(work,work))
          if(.not.ieee_is_finite(orthogonal_norm)) then
            status=coordinate
            go to 900
          end if
          hessenberg(j+1,j,coordinate)=orthogonal_norm
          if(orthogonal_norm>100.0_dp*epsilon(1.0_dp)) &
            basis(:,j+1,coordinate)=work/orthogonal_norm

          do i=1,j-1
            temp=givens_c(i,coordinate)*hessenberg(i,j,coordinate)+ &
              givens_s(i,coordinate)*hessenberg(i+1,j,coordinate)
            hessenberg(i+1,j,coordinate)= &
              -givens_s(i,coordinate)*hessenberg(i,j,coordinate)+ &
              givens_c(i,coordinate)*hessenberg(i+1,j,coordinate)
            hessenberg(i,j,coordinate)=temp
          end do
          denominator=hypot(hessenberg(j,j,coordinate), &
            hessenberg(j+1,j,coordinate))
          if(.not.ieee_is_finite(denominator)) then
            status=coordinate
            go to 900
          end if
          residual_scale=tolerance*max(1.0_dp,norm_rhs(coordinate))
          if(denominator<=100.0_dp*epsilon(1.0_dp)) then
            if(abs(g(j+1,coordinate))>residual_scale) then
              status=coordinate
              go to 900
            end if
            givens_c(j,coordinate)=1.0_dp
            givens_s(j,coordinate)=0.0_dp
          else
            givens_c(j,coordinate)=hessenberg(j,j,coordinate)/denominator
            givens_s(j,coordinate)=hessenberg(j+1,j,coordinate)/denominator
          end if
          hessenberg(j,j,coordinate)= &
            givens_c(j,coordinate)*hessenberg(j,j,coordinate)+ &
            givens_s(j,coordinate)*hessenberg(j+1,j,coordinate)
          hessenberg(j+1,j,coordinate)=0.0_dp
          temp=givens_c(j,coordinate)*g(j,coordinate)+ &
            givens_s(j,coordinate)*g(j+1,coordinate)
          g(j+1,coordinate)= &
            -givens_s(j,coordinate)*g(j,coordinate)+ &
            givens_c(j,coordinate)*g(j+1,coordinate)
          g(j,coordinate)=temp
          used(coordinate)=j
          if(abs(g(j+1,coordinate))<=residual_scale) &
            cycle_done(coordinate)=.true.
        end do
        iteration=iteration+1
      end do

      do coordinate=1,ncoord
        if(converged(coordinate)) cycle
        if(used(coordinate)<=0) then
          status=coordinate
          go to 900
        end if
        y(1:used(coordinate),coordinate)= &
          g(1:used(coordinate),coordinate)
        do i=used(coordinate),1,-1
          denominator=hessenberg(i,i,coordinate)
          if(abs(denominator)<=100.0_dp*epsilon(1.0_dp) .or. &
             .not.ieee_is_finite(denominator)) then
            status=coordinate
            go to 900
          end if
          if(i<used(coordinate)) y(i,coordinate)=y(i,coordinate)- &
            dot_product(hessenberg(i,i+1:used(coordinate),coordinate), &
              y(i+1:used(coordinate),coordinate))
          y(i,coordinate)=y(i,coordinate)/denominator
        end do
        solution(:,coordinate)=solution(:,coordinate)+ &
          matmul(preconditioned_basis(:,1:used(coordinate),coordinate), &
            y(1:used(coordinate),coordinate))
        if(any(.not.ieee_is_finite(solution(:,coordinate)))) then
          status=coordinate
          go to 900
        end if
      end do
      if(iteration==cycle_start) exit
    end do

    call apply_operator(solution,applied,operator_status)
    if(operator_status/=0 .or. any(.not.ieee_is_finite(applied))) then
      status=-3
      go to 900
    end if
    residual_max=0.0_dp
    do coordinate=1,ncoord
      residual(:,coordinate)=applied(:,coordinate)-rhs(:,coordinate)
      current_norm=sqrt(dot_product(residual(:,coordinate), &
        residual(:,coordinate)))
      component_residual=maxval(abs(residual(:,coordinate)))
      residual_max=max(residual_max,component_residual)
      if(.not.ieee_is_finite(current_norm) .or. &
         current_norm>tolerance*max(1.0_dp,norm_rhs(coordinate))) then
        status=coordinate
        go to 900
      end if
    end do

900 continue
    deallocate(basis,preconditioned_basis,hessenberg,givens_c,givens_s,g,y, &
      residual,applied,active_vectors,active_preconditioned,active_images, &
      work,norm_rhs,active_columns,used,converged,cycle_done)
  end subroutine solve_batched_gmres

!###############################################################################

  subroutine solve_shared_block_subspace(apply_operator,rhs,solution, &
      residual_max,status,tolerance,max_iterations,apply_preconditioner)
    procedure(mrsf_tda_batch_operator) :: apply_operator
    procedure(mrsf_tda_batch_operator), optional :: apply_preconditioner
    real(kind=dp), intent(in) :: rhs(:,:),tolerance
    real(kind=dp), intent(out) :: solution(:,:),residual_max
    integer, intent(out) :: status
    integer, intent(in) :: max_iterations

    integer, parameter :: maximum_stored_vectors=256
    real(kind=dp), allocatable :: residual(:,:),trial(:,:),preconditioned(:,:), &
      new_vectors(:,:),new_images(:,:),basis(:,:),images(:,:),cycle_rhs(:,:), &
      base_solution(:,:),coefficients(:,:),applied(:,:),norm_rhs(:)
    real(kind=dp) :: candidate_norm,component_residual,current_norm, &
      orthogonal_norm,threshold
    integer, allocatable :: active_columns(:)
    logical, allocatable :: converged(:)
    integer :: available,coordinate,iteration,ls_status,n,nactive,nbasis, &
      ncoord,new_count,operator_status,pass,previous_basis,trial_column

    n=size(rhs,1)
    ncoord=size(rhs,2)
    status=0
    residual_max=0.0_dp
    solution=0.0_dp
    if(n<=0 .or. ncoord<=0 .or. any(shape(solution)/=[n,ncoord]) .or. &
       .not.ieee_is_finite(tolerance) .or. tolerance<=0.0_dp .or. &
       max_iterations<=0 .or. any(.not.ieee_is_finite(rhs))) then
      status=-1
      return
    end if
    nbasis=min(n,maximum_stored_vectors)
    allocate(residual(n,ncoord),trial(n,ncoord), &
      preconditioned(n,ncoord),new_vectors(n,nbasis),new_images(n,nbasis), &
      basis(n,nbasis),images(n,nbasis),cycle_rhs(n,ncoord), &
      base_solution(n,ncoord),coefficients(nbasis,ncoord), &
      applied(n,ncoord),norm_rhs(ncoord),active_columns(ncoord), &
      converged(ncoord))
    residual=rhs
    basis=0.0_dp
    images=0.0_dp
    solution=0.0_dp
    nbasis=0
    do coordinate=1,ncoord
      norm_rhs(coordinate)=sqrt(dot_product(rhs(:,coordinate), &
        rhs(:,coordinate)))
    end do
    if(any(.not.ieee_is_finite(norm_rhs))) then
      status=-2
      go to 900
    end if
    converged=norm_rhs<=tolerance
    if(all(converged)) go to 800

    do iteration=1,max_iterations
      if(nbasis==0) then
        base_solution=solution
        cycle_rhs=residual
      end if

      call gather_active(residual,nactive)
      if(nactive==0) exit
      if(present(apply_preconditioner)) then
        call apply_preconditioner(trial(:,1:nactive), &
          preconditioned(:,1:nactive),operator_status)
        if(operator_status/=0 .or. &
           any(.not.ieee_is_finite(preconditioned(:,1:nactive)))) then
          write(error_unit,'(A,I0,A,I0)') &
            'Shared block response preconditioner failed at iteration ', &
            iteration,' with status ',operator_status
          status=-3
          go to 900
        end if
      else
        preconditioned(:,1:nactive)=trial(:,1:nactive)
      end if

      available=size(basis,2)-nbasis
      new_count=0
      do trial_column=1,nactive
        if(new_count>=available) exit
        new_vectors(:,new_count+1)=preconditioned(:,trial_column)
        candidate_norm=sqrt(dot_product(new_vectors(:,new_count+1), &
          new_vectors(:,new_count+1)))
        if(.not.ieee_is_finite(candidate_norm)) then
          status=-4
          go to 900
        end if
        do pass=1,2
          do coordinate=1,nbasis
            new_vectors(:,new_count+1)=new_vectors(:,new_count+1)- &
              basis(:,coordinate)*dot_product(basis(:,coordinate), &
                new_vectors(:,new_count+1))
          end do
          do coordinate=1,new_count
            new_vectors(:,new_count+1)=new_vectors(:,new_count+1)- &
              new_vectors(:,coordinate)*dot_product( &
                new_vectors(:,coordinate),new_vectors(:,new_count+1))
          end do
        end do
        orthogonal_norm=sqrt(dot_product(new_vectors(:,new_count+1), &
          new_vectors(:,new_count+1)))
        threshold=100.0_dp*epsilon(1.0_dp)*max(1.0_dp,candidate_norm)
        if(orthogonal_norm<=threshold) cycle
        new_count=new_count+1
        new_vectors(:,new_count)=new_vectors(:,new_count)/orthogonal_norm
      end do
      if(new_count==0) then
        write(error_unit,'(A,I0,A,I0,A,I0)') &
          'Shared block response generated no independent vector at iteration ', &
          iteration,', stored basis ',nbasis,', active RHS ',nactive
        status=first_unconverged()
        go to 900
      end if

      call apply_operator(new_vectors(:,1:new_count), &
        new_images(:,1:new_count),operator_status)
      if(operator_status/=0 .or. &
         any(.not.ieee_is_finite(new_images(:,1:new_count)))) then
        write(error_unit,'(A,I0,A,I0,A,I0)') &
          'Shared block response operator failed at iteration ',iteration, &
          ' with status ',operator_status,', block width ',new_count
        status=-3
        go to 900
      end if
      previous_basis=nbasis
      basis(:,nbasis+1:nbasis+new_count)=new_vectors(:,1:new_count)
      images(:,nbasis+1:nbasis+new_count)=new_images(:,1:new_count)
      nbasis=nbasis+new_count

      call solve_block_least_squares(images(:,1:nbasis),cycle_rhs, &
        coefficients(1:nbasis,:),ls_status)
      if(ls_status/=0) then
        write(error_unit,'(A,I0,A,I0)') &
          'Shared response subspace least-squares failure at basis ', &
          previous_basis+ls_status,' of ',nbasis
        status=-5
        go to 900
      end if
      solution=base_solution+matmul(basis(:,1:nbasis), &
        coefficients(1:nbasis,:))
      residual=cycle_rhs-matmul(images(:,1:nbasis), &
        coefficients(1:nbasis,:))
      do coordinate=1,ncoord
        current_norm=sqrt(dot_product(residual(:,coordinate), &
          residual(:,coordinate)))
        if(.not.ieee_is_finite(current_norm)) then
          status=coordinate
          go to 900
        end if
        converged(coordinate)=current_norm<= &
          tolerance*max(1.0_dp,norm_rhs(coordinate))
      end do
      if(all(converged)) exit
      if(nbasis==size(basis,2)) nbasis=0
    end do
    if(.not.all(converged)) then
      do coordinate=1,ncoord
        if(.not.converged(coordinate)) then
          write(error_unit,'(A,I0,A,I0)') &
            'Shared block response failed for coordinate ',coordinate, &
            ' after block iterations ',max_iterations
          status=coordinate
          exit
        end if
      end do
      go to 900
    end if

800 continue
    if(all(norm_rhs<=tolerance)) then
      residual_max=0.0_dp
      go to 900
    end if
    call apply_operator(solution,applied,operator_status)
    if(operator_status/=0 .or. any(.not.ieee_is_finite(applied))) then
      write(error_unit,'(A,I0)') &
        'Shared block response final operator certification failed; status ', &
        operator_status
      status=-3
      go to 900
    end if
    do coordinate=1,ncoord
      current_norm=sqrt(dot_product(applied(:,coordinate)-rhs(:,coordinate), &
        applied(:,coordinate)-rhs(:,coordinate)))
      component_residual=maxval(abs(applied(:,coordinate)-rhs(:,coordinate)))
      residual_max=max(residual_max,component_residual)
      if(.not.ieee_is_finite(current_norm) .or. &
         component_residual>max(1.0e-8_dp,100.0_dp*tolerance)* &
           max(1.0_dp,norm_rhs(coordinate))) then
        write(error_unit,'(A,I0,A,1P,E12.4,A,E12.4,A,E12.4)') &
          'Shared block response final residual failed for RHS ',coordinate, &
          ': max=',component_residual,', norm=',current_norm, &
          ', reference norm=',norm_rhs(coordinate)
        status=coordinate
        go to 900
      end if
    end do

900 continue
    deallocate(residual,trial,preconditioned,new_vectors,new_images,basis, &
      images,cycle_rhs,base_solution,coefficients,applied,norm_rhs, &
      active_columns,converged)

  contains

    subroutine gather_active(vectors,nfound)
      real(kind=dp), intent(in) :: vectors(:,:)
      integer, intent(out) :: nfound
      integer :: local_coordinate
      nfound=0
      do local_coordinate=1,ncoord
        if(.not.converged(local_coordinate)) then
          nfound=nfound+1
          active_columns(nfound)=local_coordinate
          trial(:,nfound)=vectors(:,local_coordinate)
        end if
      end do
    end subroutine gather_active

    integer function first_unconverged()
      integer :: local_coordinate
      first_unconverged=1
      do local_coordinate=1,ncoord
        if(.not.converged(local_coordinate)) then
          first_unconverged=local_coordinate
          return
        end if
      end do
    end function first_unconverged

    subroutine solve_block_least_squares(operator_images,block_rhs, &
                                         block_coefficients,local_status)
      real(kind=dp), intent(in) :: operator_images(:,:),block_rhs(:,:)
      real(kind=dp), intent(out) :: block_coefficients(:,:)
      integer, intent(out) :: local_status
      real(kind=dp), allocatable :: q(:,:),r(:,:),reduced_rhs(:,:)
      real(kind=dp) :: diagonal,overlap,scale
      integer :: column,local_pass,previous,rhs_column,subspace_size

      subspace_size=size(operator_images,2)
      local_status=0
      block_coefficients=0.0_dp
      if(size(operator_images,1)/=n .or. &
         any(shape(block_rhs)/=[n,ncoord]) .or. &
         any(shape(block_coefficients)/=[subspace_size,ncoord])) then
        local_status=-1
        return
      end if
      allocate(q(n,subspace_size),r(subspace_size,subspace_size), &
        reduced_rhs(subspace_size,ncoord),source=0.0_dp)
      q=operator_images
      do column=1,subspace_size
        scale=sqrt(dot_product(q(:,column),q(:,column)))
        do local_pass=1,2
          do previous=1,column-1
            overlap=dot_product(q(:,previous),q(:,column))
            r(previous,column)=r(previous,column)+overlap
            q(:,column)=q(:,column)-overlap*q(:,previous)
          end do
        end do
        diagonal=sqrt(dot_product(q(:,column),q(:,column)))
        if(.not.ieee_is_finite(diagonal) .or. &
           diagonal<=100.0_dp*epsilon(1.0_dp)*max(1.0_dp,scale)) then
          local_status=column
          deallocate(q,r,reduced_rhs)
          return
        end if
        r(column,column)=diagonal
        q(:,column)=q(:,column)/diagonal
      end do
      reduced_rhs=matmul(transpose(q),block_rhs)
      do rhs_column=1,ncoord
        do column=subspace_size,1,-1
          block_coefficients(column,rhs_column)= &
            (reduced_rhs(column,rhs_column)-dot_product( &
              r(column,column+1:subspace_size), &
              block_coefficients(column+1:subspace_size,rhs_column)))/ &
            r(column,column)
        end do
      end do
      if(any(.not.ieee_is_finite(block_coefficients))) local_status=-2
      deallocate(q,r,reduced_rhs)
    end subroutine solve_block_least_squares

  end subroutine solve_shared_block_subspace

!###############################################################################

  subroutine apply_coupled_matrix(amb, apb, omega, x, y, nexc)
    real(kind=dp), intent(in) :: amb(:,:), apb(:,:), omega, x(:)
    real(kind=dp), intent(out) :: y(:)
    integer, intent(in) :: nexc

    y(1:nexc) = matmul(amb, x(1:nexc)) - omega*x(nexc+1:2*nexc)
    y(nexc+1:2*nexc) = matmul(apb, x(nexc+1:2*nexc)) - &
                       omega*x(1:nexc)
  end subroutine apply_coupled_matrix

!###############################################################################

  subroutine project_from_null(vector, nullvec, norm2)
    real(kind=dp), intent(inout) :: vector(:)
    real(kind=dp), intent(in) :: nullvec(:), norm2

    vector = vector - nullvec*(dot_product(nullvec, vector)/norm2)
  end subroutine project_from_null

end module tdhf_hessian_response_mod
