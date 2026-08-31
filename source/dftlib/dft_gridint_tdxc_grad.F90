module mod_dft_gridint_tdxc_grad

  use precision, only: fp
  use mod_dft_gridint, only: xc_engine_t, xc_consumer_t
  use mod_dft_gridint, only: OQP_FUNTYP_LDA, OQP_FUNTYP_GGA, OQP_FUNTYP_MGGA
  use mod_dft_gridint, only: compAtGradRho, compAtGradDRho, compAtGradTau

  implicit none

!-------------------------------------------------------------------------------

  type, extends(xc_consumer_t) :: xc_consumer_tdg_t
    integer :: nMtx = 1
    logical :: do_fxc = .true. !< Whether to compute dF_xc / dR_i
    logical :: do_ground_state = .true. !< Whether to add g.s. XC gradient contribution
    logical :: do_weight_derivative = .false. !< Moving-grid response for a linear probe
    integer :: part_fun_type = 0
    logical :: has_surface_shift = .false.
    real(kind=fp), allocatable :: atom_xyz(:,:)
    logical, allocatable :: dummy_atom(:)
    real(kind=fp), allocatable :: surface_shift(:,:)
    real(kind=fp), allocatable :: part_rij(:,:)
    real(kind=fp), allocatable :: part_rhat(:,:,:)
    real(kind=fp), allocatable :: part_dist(:,:)
    real(kind=fp), allocatable :: part_cells(:,:)
    real(kind=fp), allocatable :: part_dlog(:,:,:,:)
    real(kind=fp), allocatable :: part_dsum(:,:,:)
    real(kind=fp), pointer :: pa(:,:,:)
    real(kind=fp), pointer :: pb(:,:,:)
    real(kind=fp), pointer :: xa(:,:,:)
    real(kind=fp), pointer :: xb(:,:,:)
    real(kind=fp), allocatable :: rrho(:,:,:,:)
    real(kind=fp), allocatable :: drrho(:,:,:,:,:)
    real(kind=fp), allocatable :: rtau(:,:,:,:)
    real(kind=fp), allocatable :: bfgrad(:,:,:)
    real(kind=fp), allocatable :: nucgrad(:,:,:)
    real(kind=fp), allocatable :: probe_value(:,:)
    real(kind=fp), allocatable :: grad_d(:,:,:,:) !< density gradient
    real(kind=fp), allocatable :: grad_p(:,:,:,:) !< diff. density gradient
    real(kind=fp), allocatable :: grad_x(:,:,:,:) !< transition (X+Y) gradient
!   Temporary storage
    real(kind=fp), allocatable :: tmpGrad_(:,:)
    real(kind=fp), allocatable :: tmp_(:,:,:,:)
    real(kind=fp), allocatable :: tmpV_(:,:,:)
    real(kind=fp), allocatable :: tmpG1_(:,:,:)

  contains
    procedure :: parallel_start
    procedure :: parallel_stop
    procedure :: resetGradPointers
    procedure :: resetPointers
    procedure :: update
    procedure :: postUpdate
    procedure :: clean
  end type

!-------------------------------------------------------------------------------

  private
  public tddft_xc_gradient
  public utddft_xc_gradient

!-------------------------------------------------------------------------------

contains

!-------------------------------------------------------------------------------

  subroutine parallel_start(self, xce, nthreads)
    implicit none
    class(xc_consumer_tdg_t), target, intent(inout) :: self
    class(xc_engine_t), intent(in) :: xce
    integer, intent(in) :: nthreads
    integer :: nspin, nterms, nDeriv, nat, i, j
    call clean_work(self)
    nterms = 1
    if (xce%funTyp /= OQP_FUNTYP_LDA) nterms = nterms + 3
    if (xce%funTyp == OQP_FUNTYP_MGGA) nterms = nterms + 1

    nspin = merge(2, 1, xce%hasBeta)
    nDeriv = merge(2, 1, self%do_fxc)
    allocate( &
        self%bfgrad(xce%numAOs, 3, nthreads) &
      , self%nucgrad(3, xce%numAtoms, nthreads) &
      , self%probe_value(xce%maxPts, nthreads) &
      , self%rrho(nspin, xce%maxPts, self%nMtx, nthreads) &
      , self%drrho(3, nspin, xce%maxPts, self%nMtx, nthreads) &
      , self%grad_d(xce%maxPts, nterms, nspin, nthreads) &
      , self%grad_p(xce%maxPts, nterms, nspin, nthreads) &
!   Temporary storage
      , self%tmpGrad_(xce%numAOs*3, nthreads) &
      , self%tmp_(xce%numAOs * xce%numAOs * self%nMtx, nspin, nDeriv, nthreads) &
      , self%tmpV_(xce%numAOs * xce%maxPts * self%nMtx * nspin, nDeriv, nthreads) &
      , self%tmpG1_(xce%numAOs * xce%maxPts * 3 * self%nMtx * nspin, nDeriv, nthreads) &
      , source=0.0d0)

    if (self%do_fxc) then
      allocate( &
          self%grad_x(xce%maxPts, nterms, nspin, nthreads) &
        , source=0.0d0)
    end if

    if (xce%funTyp == OQP_FUNTYP_MGGA) then
        allocate( &
            self%rtau(nSpin, xce%maxPts, self%nMtx, nthreads) &
          , source=0.0d0)
    end if

    if (self%do_weight_derivative) then
      nat = xce%numAtoms
      allocate( &
          self%part_rij(nat,nat) &
        , self%part_rhat(3,nat,nat) &
        , self%part_dist(nat,nthreads) &
        , self%part_cells(nat,nthreads) &
        , self%part_dlog(3,nat,nat,nthreads) &
        , self%part_dsum(3,nat,nthreads) &
        , source=0.0_fp)
      do i = 1, nat
        do j = 1, nat
          if (i == j) cycle
          self%part_rij(i,j) = norm2( &
            self%atom_xyz(:,i)-self%atom_xyz(:,j))
          if (self%part_rij(i,j) > tiny(1.0_fp)) &
            self%part_rhat(:,i,j) = &
              (self%atom_xyz(:,i)-self%atom_xyz(:,j))/self%part_rij(i,j)
        end do
      end do
    end if
  end subroutine

!-------------------------------------------------------------------------------

  subroutine parallel_stop(self)
    implicit none
    class(xc_consumer_tdg_t), intent(inout) :: self

    if (ubound(self%bfGrad,3) /= 1) then
      self%bfGrad(:,:,lbound(self%bfGrad,3)) = sum(self%bfGrad, dim=3)
      self%nucGrad(:,:,lbound(self%nucGrad,3)) = sum(self%nucGrad, dim=3)
    end if
    call self%pe%allreduce(self%bfGrad(:,:,1), &
              size(self%bfGrad(:,:,1)))
    call self%pe%allreduce(self%nucGrad(:,:,1), &
              size(self%nucGrad(:,:,1)))
  end subroutine

!-------------------------------------------------------------------------------

  subroutine clean_work(self)
    implicit none
    class(xc_consumer_tdg_t), intent(inout) :: self
    if (allocated(self%bfgrad)) deallocate(self%bfgrad)
    if (allocated(self%nucgrad)) deallocate(self%nucgrad)
    if (allocated(self%probe_value)) deallocate(self%probe_value)
    if (allocated(self%rrho)) deallocate(self%rrho)
    if (allocated(self%drrho)) deallocate(self%drrho)
    if (allocated(self%rtau)) deallocate(self%rtau)
    if (allocated(self%grad_d)) deallocate(self%grad_d)
    if (allocated(self%grad_p)) deallocate(self%grad_p)
    if (allocated(self%grad_x)) deallocate(self%grad_x)
    if (allocated(self%tmpGrad_)) deallocate(self%tmpGrad_)
    if (allocated(self%tmp_)) deallocate(self%tmp_)
    if (allocated(self%tmpV_)) deallocate(self%tmpV_)
    if (allocated(self%tmpG1_)) deallocate(self%tmpG1_)
    if (allocated(self%part_rij)) deallocate(self%part_rij)
    if (allocated(self%part_rhat)) deallocate(self%part_rhat)
    if (allocated(self%part_dist)) deallocate(self%part_dist)
    if (allocated(self%part_cells)) deallocate(self%part_cells)
    if (allocated(self%part_dlog)) deallocate(self%part_dlog)
    if (allocated(self%part_dsum)) deallocate(self%part_dsum)
  end subroutine clean_work

!-------------------------------------------------------------------------------

  subroutine clean(self)
    implicit none
    class(xc_consumer_tdg_t), intent(inout) :: self
    call clean_work(self)
    if (allocated(self%atom_xyz)) deallocate(self%atom_xyz)
    if (allocated(self%dummy_atom)) deallocate(self%dummy_atom)
    if (allocated(self%surface_shift)) deallocate(self%surface_shift)
  end subroutine

!-------------------------------------------------------------------------------
!> @brief Adjust internal memory storage for a given
!>  number of pruned grid points
!> @author Konstantin Komarov
 subroutine resetGradPointers(self, xce, tmpGrad, tmpV, tmpG1, myThread)
    class(xc_consumer_tdg_t), target, intent(inout) :: self
    class(xc_engine_t), intent(in) :: xce
    real(kind=fp), intent(out), pointer :: tmpGrad(:,:)
    real(kind=fp), intent(out), pointer, optional :: tmpV(:,:,:,:,:)
    real(kind=fp), intent(out), pointer, optional :: tmpG1(:,:,:,:,:,:)
    integer, intent(in) :: myThread
    integer :: nSpin

    nspin = merge(2, 1, xce%hasBeta)
    associate ( numAOs => xce%numAOs_p &  ! number of pruned AOs
              , numPts => xce%numPts &
              , nMtx   => self%nMtx &
      )

      tmpGrad(1:numAOs,1:3) => self%tmpGrad_(1:numAOs*3, myThread)

      if (present(tmpV)) &
        tmpV(1:numAOs, 1:numPts, 1:nMtx, 1:nSpin, 1:1) => &
           self%tmpV_(1:numAOs*numPts*nMtx*nspin, 1, myThread)

      if (present(tmpG1)) &
        tmpG1(1:numAOs, 1:numPts, 1:3, 1:nMtx, 1:nSpin, 1:1) => &
          self%tmpG1_(1:numAOs*numPts*3*nMtx*nspin, 1, mythread)

      if (present(tmpV) .and. self%do_fxc) &
        tmpV(1:numAOs, 1:numPts, 1:nMtx, 1:nSpin, 2:2) => &
          self%tmpV_(1:numAOs*numPts*nMtx*nspin, 2, myThread)

      if (present(tmpG1) .and. self%do_fxc) &
          tmpG1(1:numAOs, 1:numPts, 1:3, 1:nMtx, 1:nSpin, 2:2) => &
            self%tmpG1_(1:numAOs*numPts*3*nMtx*nspin, 2, mythread)

    end associate

 end subroutine

 subroutine resetPointers(self, xce, Pa, Pb, Xa, Xb, &
         Pa_p, Pb_p, Xa_p, Xb_p,  myThread)
    class(xc_consumer_tdg_t), target, intent(inout) :: self
    class(xc_engine_t), intent(in) :: xce
    real(kind=fp), intent(in), target :: Pa(:,:,:)
    real(kind=fp), intent(in), target :: Pb(:,:,:)
    real(kind=fp), intent(in), target :: Xa(:,:,:)
    real(kind=fp), intent(in), target :: Xb(:,:,:)
    real(kind=fp), intent(out), pointer :: Pa_p(:,:,:)  ! pruned
    real(kind=fp), intent(out), pointer :: Pb_p(:,:,:)  ! pruned
    real(kind=fp), intent(out), pointer :: Xa_p(:,:,:)  ! pruned
    real(kind=fp), intent(out), pointer :: Xb_p(:,:,:)  ! pruned
    integer, intent(in) :: myThread
    integer :: nSpin

    nspin = merge(2, 1, xce%hasBeta)
    associate ( indices => xce%indices_p &
              , numAOs  => xce%numAOs_p &  ! number of pruned numAOs
              , numPts  => xce%numPts &
              , nMtx    => self%nMtx &
      )
      if (xce%skip_p) then
!       no pruned AOs
        Pa_p => Pa
        if (xce%hasBeta) Pb_p => Pb
        if (self%do_fxc) Xa_p => Xa
        if (self%do_fxc .and. xce%hasBeta) Xb_p => Xb

      else

!       pruned AOs
        Pa_p(1:numAOs, 1:numAOs, 1:nMtx) => self%tmp_(1:numAOs*numAOs*nMtx, 1, 1, myThread)
!       Compress matrix
        Pa_p(1:numAOs, 1:numAOs,:) = Pa(indices(1:numAOs), indices(1:numAOs),:)

        ! if do dF_xc / dR_i
        if (self%do_fxc) then
          Xa_p(1:numAOs, 1:numAOs, 1:nMtx) => self%tmp_(1:numAOs*numAOs*nMtx, 1, 2, myThread)
          Xa_p(1:numAOs, 1:numAOs,:) = Xa(indices(1:numAOs), indices(1:numAOs),:)
        end if

        if (xce%hasBeta) then
          Pb_p(1:numAOs, 1:numAOs, 1:nMtx) => self%tmp_(1:numAOs*numAOs*nMtx, 2, 1, myThread)
          Pb_p(1:numAOs, 1:numAOs,:) = Pb(indices(1:numAOs), indices(1:numAOs),:)
          if (self%do_fxc) then
            Xb_p(1:numAOs, 1:numAOs, 1:nMtx) => self%tmp_(1:numAOs*numAOs*nMtx, 2, 2, myThread)
            Xb_p(1:numAOs, 1:numAOs,:) = Xb(indices(1:numAOs), indices(1:numAOs),:)
          end if
        end if

      end if

    end associate

 end subroutine

 subroutine update(self, xce, mythread)

    class(xc_consumer_tdg_t), intent(inout) :: self
    class(xc_engine_t), intent(in) :: xce
    integer :: mythread
    real(kind=fp), pointer :: tmpGrad(:,:)
    real(kind=fp), pointer :: Pa(:,:,:)
    real(kind=fp), pointer :: Pb(:,:,:)
    real(kind=fp), pointer :: Xa(:,:,:)
    real(kind=fp), pointer :: Xb(:,:,:)
    real(kind=fp), pointer :: tmpV(:,:,:,:,:)
    real(kind=fp), pointer :: tmpG1(:,:,:,:,:,:)

    call self%resetGradPointers(xce, tmpGrad, tmpV, tmpG1, myThread)

    ! Needs to nullify it for each update
    tmpGrad = 0.0d0

    associate ( bfgrad => self%bfgrad(:,:,mythread) &
              , grad_d => self%grad_d(:,:,:,mythread) &
              , grad_p => self%grad_p(:,:,:,mythread) &
              , grad_x => self%grad_x(:,:,:,mythread) &
              , aoV    => xce%aoV &
              , aoG1   => xce%aoG1 &
              , aoG2   => xce%aoG2 &
              , moVA   => xce%moVA &
              , moVB   => xce%moVB &
              , moG1A  => xce%moG1A &
              , moG1B  => xce%moG1B &
              , rrho   => self%rrho(:,:,:,mythread)  &
              , drrho  => self%drrho(:,:,:,:,mythread)  &
              , rtau   => self%rtau(:,:,:,mythread)  &
              , drho   => xce%xclib%drho  &
              , numPts => xce%numPts &
              , xc     => xce%XCLib &
              , ids    => xce%XCLib%ids &
      )

      ! Compute "MOs" which correspond to the ground
      ! state and the relaxed difference density matrices
      ! They are basically right sides of
      ! \Phi (A \Phi), and \Phi (A \nabla\Psi)
      ! where A is some density-like matrix
      call self%resetPointers(xce, self%pa, self%pb, self%xa, self%xb, &
                Pa, Pb, Xa, Xb, myThread)

      call xce%compRMOs(Pa, tmpV(:,:,:,1,1))
      call xce%compRMOGs(Pa, tmpG1(:,:,:,:,1,1))

      if (xce%hasBeta) then
        call xce%compRMOs(Pb, tmpV(:,:,:,2,1))
        call xce%compRMOGs(Pb, tmpG1(:,:,:,:,2,1))
      end if

      ! d V_xc / d R_i
      ! Compute difference densities: \rho, \nabla\rho, and \tau
      call xce%compRRho(tmpV(:,:,:,:,1), rRho)
      call xce%compRDRho(tmpV(:,:,:,:,1), drRho)
      if (xce%funTyp == OQP_FUNTYP_MGGA) then
        call xce%compRTau(tmpG1(:,:,:,:,:,1), rTau)
      end if
      ! Compute XC 1st and 2nd derivative and compute terms for contraction
      ! with g.s. density (`grad_d`) and difference density (`grad_p`)
      if (xce%hasBeta) then
        call grad_v_xc(self, xce, mythread)
      else
        call grad_v_xc_np(self, xce, mythread)
      end if

      ! d F_xc / dR_i
      if (self%do_fxc) then
        ! Compute "MOs" which correspond to the `X+Y` transition density
        call xce%compRMOs(Xa, tmpV(:,:,:,1,2))
        call xce%compRMOGs(Xa, tmpG1(:,:,:,:,1,2))
        if (xce%hasBeta) then
          call xce%compRMOs(Xb, tmpV(:,:,:,2,2))
          call xce%compRMOGs(Xb, tmpG1(:,:,:,:,2,2))
        end if
        ! Compute transition densities: \rho, \nabla\rho, and \tau
        call xce%compRRho(tmpV(:,:,:,:,2), rrho)
        call xce%compRDRho(tmpV(:,:,:,:,2), drrho)
        if (xce%funTyp == OQP_FUNTYP_MGGA) then
          call xce%compRTau(tmpG1(:,:,:,:,:,2), rtau)
        end if
        ! Compute XC 1st-3rd derivatives and compute terms for contraction
        ! with g.s. density (`grad_d`) and transition density (`grad_x`)
        if (xce%hasBeta) then
          call grad_f_xc(self, xce, mythread)
        else
          call grad_f_xc_np(self, xce, mythread)
        end if
      end if

      if (self%do_ground_state) then
        ! `grad_p` also contains functional derivative terms
        ! required to compute ground state contribution to the gradient.
        ! We can do it here instead of separate G.S. DFT gradient run
        ! To contract them with the ground state density
        ! we add it here to the `grad_d`
        grad_d = grad_d + grad_p
      end if
      ! In RHF case g.s. density is twice the alpha density
      ! TODO: make it consistent
      if (.not. xce%hasBeta) grad_d = 0.5*grad_d

      ! Compute contribution to the AO gradient from all `grad_P` terms
      call compAtGradAll(tmpGrad, grad_D(:,:,1), xce%funTyp, &
                         moVA, moG1A, aoG1, aoG2, numPts)
      call compAtGradAll(tmpGrad, grad_P(:,:,1), xce%funTyp, &
                         tmpV(:,:,1,1,1), tmpG1(:,:,:,1,1,1), aoG1, aoG2, numPts)

      if (xce%hasBeta) then
        call compAtGradAll(tmpGrad, grad_D(:,:,2), xce%funTyp, &
                           moVB, moG1B, aoG1, aoG2, numPts)
        call compAtGradAll(tmpGrad, grad_P(:,:,2), xce%funTyp, &
                           tmpV(:,:,1,2,1), tmpG1(:,:,:,1,2,1), aoG1, aoG2, numPts)
      end if

      ! d F_xc / dR_i
      if (self%do_fxc) then

        ! Compute contribution to the AO gradient from all `grad_X` terms
        call compAtGradAll(tmpGrad, grad_X(:,:,1), xce%funTyp, &
                           tmpV(:,:,1,1,2), tmpG1(:,:,:,1,1,2), aoG1, aoG2, numPts)
        if (xce%hasBeta) &
          call compAtGradAll(tmpGrad, grad_X(:,:,2), xce%funTyp, &
                             tmpV(:,:,1,2,2), tmpG1(:,:,:,1,2,2), aoG1, aoG2, numPts)

      end if

      if (self%do_weight_derivative) then
        call add_partition_weight_gradient(self, xce, mythread)
        ! Differentiate the discrete atom-centred quadrature consistently:
        ! the grid point moves with the atom that owns the current slice.
        self%nucgrad(:,xce%currAtom,mythread) = &
          self%nucgrad(:,xce%currAtom,mythread) + sum(tmpGrad, dim=1)
      end if

   end associate

 end subroutine

!> Add the derivative of normalized atom-centred fuzzy-cell weights for a
!> linear density probe P.  The AO/basis contribution is accumulated in
!> tmpGrad; this routine supplies q_P d(log p_owner)/dR.  The owner-motion
!> contribution is added by update immediately after this call.
 subroutine add_partition_weight_gradient(self, xce, mythread)
    use mod_dft_partfunc, only: partition_function

    class(xc_consumer_tdg_t), intent(inout) :: self
    class(xc_engine_t), intent(in) :: xce
    integer, intent(in) :: mythread

    type(partition_function) :: partfunc
    real(kind=fp) :: point(3), ui(3), uj(3), df(3)
    real(kind=fp) :: dri(3), drj(3), drij(3), dmu(3)
    real(kind=fp) :: mu0, mu, f, fi, fj, dfi_scale, dfj_scale
    real(kind=fp) :: sumc, p_owner, q_weighted, aij, dfactor
    integer :: nat, owner, ipt, i, j, b, ib, nb
    integer :: derivative_atoms(3)

    nat = xce%numAtoms
    owner = xce%currAtom
    if (owner < 1 .or. owner > nat) return
    if (.not. allocated(self%atom_xyz)) return
    call partfunc%set(self%part_fun_type)

    associate ( &
        dist => self%part_dist(:,mythread) &
      , cells => self%part_cells(:,mythread) &
      , dlog => self%part_dlog(:,:,:,mythread) &
      , dsum => self%part_dsum(:,:,mythread))
      do ipt = 1, xce%numPts
        q_weighted = self%probe_value(ipt,mythread)
        if (abs(q_weighted) <= tiny(1.0_fp)) cycle

        point = xce%xyzw(ipt,1:3)
        do i = 1, nat
          dist(i) = norm2(point-self%atom_xyz(:,i))
        end do

        cells = 1.0_fp
        where (self%dummy_atom) cells = 0.0_fp
        dlog = 0.0_fp

        do i = 2, nat
          if (self%dummy_atom(i)) cycle
          do j = 1, i-1
            if (self%dummy_atom(j)) cycle
            if (self%part_rij(j,i) <= tiny(1.0_fp)) cycle

            mu0 = (dist(i)-dist(j))/self%part_rij(j,i)
            mu = mu0
            aij = 0.0_fp
            if (self%has_surface_shift) then
              aij = self%surface_shift(j,i)
              mu = mu0 + aij*(1.0_fp-mu0*mu0)
            end if
            f = partfunc%eval(mu)
            fi = abs(f)
            fj = abs(1.0_fp-f)
            dfi_scale = sign(1.0_fp, f)
            dfj_scale = -sign(1.0_fp, 1.0_fp-f)

            ui = 0.0_fp
            uj = 0.0_fp
            if (dist(i) > tiny(1.0_fp)) &
              ui = (point-self%atom_xyz(:,i))/dist(i)
            if (dist(j) > tiny(1.0_fp)) &
              uj = (point-self%atom_xyz(:,j))/dist(j)

            ! Only the point owner and the two atoms in this pair affect
            ! mu_ij.  Logarithmic derivatives reduce the work to
            ! O(Ngrid*Natom**2) and avoid allocation in this loop.
            derivative_atoms = owner
            derivative_atoms(2) = i
            derivative_atoms(3) = j
            nb = 1
            if (i /= owner) nb = nb + 1
            if (j /= owner .and. j /= i) nb = nb + 1
            if (nb == 2) then
              if (i == owner) derivative_atoms(2) = j
            else if (nb == 3) then
              derivative_atoms(2) = i
              derivative_atoms(3) = j
            end if

            do ib = 1, nb
              b = derivative_atoms(ib)
              dri = ui * real(merge(1,0,b == owner) - &
                              merge(1,0,b == i), fp)
              drj = uj * real(merge(1,0,b == owner) - &
                              merge(1,0,b == j), fp)
              drij = self%part_rhat(:,i,j) * &
                      real(merge(1,0,b == i) - merge(1,0,b == j), fp)
              dmu = (dri-drj-mu0*drij)/self%part_rij(j,i)
              if (self%has_surface_shift) &
                dmu = (1.0_fp-2.0_fp*aij*mu0)*dmu
              df = partfunc%deriv(mu)*dmu
              if (fi > tiny(1.0_fp)) &
                dlog(:,b,i) = dlog(:,b,i) + dfi_scale*df/fi
              if (fj > tiny(1.0_fp)) &
                dlog(:,b,j) = dlog(:,b,j) + dfj_scale*df/fj
            end do
            cells(i) = cells(i)*fi
            cells(j) = cells(j)*fj
          end do
        end do

        sumc = sum(cells)
        if (sumc <= tiny(1.0_fp)) cycle
        p_owner = cells(owner)/sumc
        if (p_owner <= sqrt(tiny(1.0_fp))) cycle
        dsum = 0.0_fp
        do i = 1, nat
          do b = 1, nat
            dsum(:,b) = dsum(:,b) + cells(i)*dlog(:,b,i)
          end do
        end do
        do b = 1, nat
          dfactor = 1.0_fp/sumc
          self%nucgrad(:,b,mythread) = self%nucgrad(:,b,mythread) &
            + q_weighted * (dlog(:,b,owner)-dfactor*dsum(:,b))
        end do
      end do
    end associate
 end subroutine add_partition_weight_gradient

 subroutine postUpdate(self, xce, mythread)

    class(xc_consumer_tdg_t), intent(inout) :: self
    class(xc_engine_t), intent(in) :: xce
    integer :: mythread

    real(kind=fp), pointer :: tmpGrad(:,:)

    call self%resetGradPointers(xce, tmpGrad,  myThread=myThread)

    associate ( numAOs  => xce%numAOs_p &  ! number of pruned AOs
              , indices => xce%indices_p &
      )

      if (xce%skip_p) then

        self%bfGrad(:,:,myThread) = self%bfGrad(:,:,myThread) + tmpGrad

      else

        self%bfGrad(indices(1:numAOs), :, mythread) = &
          self%bfGrad(indices(1:numAOs), :, mythread) + tmpGrad(1:numAOs, :)

      end if

   end associate
 end subroutine

!> @brief Compute contribution to the AO gradient from
!>   LDA, GGA and metaGGA functional derivatives
!> @param[inout] bfGrad   array of gradient contributions per AO
!> @param[in]    fgrad    XC gradient
!> @param[in]    funTyp   type of the XC functional
!> @param[in]    moV      MO-like orbital values
!> @param[in]    moG1     MO-like orbital gradients
!> @param[in]    aoG1     AO orbital gradients
!> @param[in]    aoG2     AO orbital 2nd derivatives
!> @param[in]    npts     number of grid points in a chunk
!> @author Vladimir Mironov
 subroutine compAtGradAll(bfGrad, fgrad, funTyp, moV, moG1, aoG1, aoG2, npts)
    integer, intent(in) :: funTyp
    real(kind=fp), intent(in) :: fgrad(:,:)
    real(kind=fp), intent(inout) :: bfGrad(:,:)
    real(kind=fp), contiguous, intent(in) :: moV(:,:), aoG1(:,:,:)
    real(kind=fp), contiguous, intent(in) :: moG1(:,:,:)
    real(kind=fp), contiguous, intent(in) :: aoG2(:,:,:)
    integer, intent(in) :: npts

!   LDA gradient
    call compAtGradRho(bfGrad, fgrad(:,1), moV(:,:), aoG1, npts)

!   GGA gradient
    if (funTyp /= OQP_FUNTYP_LDA) then
        call compAtGradDRho(bfGrad, fgrad(:,2:4), &
                moV, moG1, aoG1, aoG2, npts)
    end if

    if (funTyp == OQP_FUNTYP_MGGA) then
        call compAtGradTau(bfGrad, fgrad(:,5), &
                moG1(:,:,:), aoG2, npts)
    end if

 end subroutine

!> @brief Compute derivative terms of \sum_ij V^xc_ij P_ij
!>  w.r.t. atomic coordinates
!> @detail this subroutine update terms which should be contracted with
!>  ground state (`grad_d`) and relaxed difference densities (`grad_p`)
!> @note For further details see:
!>  [1] F.Furche, R.Ahlrichs, J.Chem.Phys. 117, p. 7433 (2002)
!>  [2] G.Scalmani et al., J.Chem.Phys. 124, 094107 (2006)
!> @author Vladimir Mironov
 subroutine grad_v_xc_np(dat, xce, mythread)
    use mod_dft_gridint, only: xc_der1, xc_der2_contr
    class(xc_engine_t) :: xce
    type(xc_consumer_tdg_t) :: dat
    integer :: mythread

    integer :: i, j
    real(kind=fp) :: d_r(2), d_s(3), d_t(2)
    real(kind=fp) :: f_r(2), f_s(3), f_t(2)
    real(kind=fp) :: rhoab(2), tauab(2), sigma(3)

    associate ( grad_d  => dat%grad_d(:,:,:,mythread) &
              , grad_p  => dat%grad_p(:,:,:,mythread) &
              , aoG1    => xce%aoG1 &
              , aoV     => xce%aoV  &
              , rrho    => dat%rrho(:,:,:,mythread)  &
              , drrho   => dat%drrho(:,:,:,:,mythread)  &
              , rtau    => dat%rtau(:,:,:,mythread)  &
              , drho    => xce%xclib%drho  &
              , numAOs  => xce%numAOs &
              , numPts  => xce%numPts &
              , nMtx    => dat%nMtx &
              , xc      => xce%XCLib &
              , ids     => xce%XCLib%ids &
              )
    do j = 1, nMtx

      do i = 1, numPts

        rhoab = rrho(1,i,j)
        sigma = 0
        tauab = 0
        if (xce%funTyp /= OQP_FUNTYP_LDA) then
          ! Compute sigma terms:
          ! \nabla\rho(D) \dot \nabla\rho(P)
          sigma = 2*dot_product(drrho(:,1,i,j), drho(1:3,i))
        end if
        if (xce%funTyp == OQP_FUNTYP_MGGA) tauab = rtau(1,i,j)

        call xc_der1(xce, .false., i, d_r, d_s, d_t)
        call xc_der2_contr(xce, .false., i, &
                rhoab, sigma, tauab, &
                f_r, f_s, f_t)

        if (dat%do_weight_derivative) then
          dat%probe_value(i,mythread) = dot_product(d_r, rhoab)
          if (xce%funTyp /= OQP_FUNTYP_LDA) &
            dat%probe_value(i,mythread) = dat%probe_value(i,mythread) &
                                        + dot_product(d_s, sigma)
          if (xce%funTyp == OQP_FUNTYP_MGGA) &
            dat%probe_value(i,mythread) = dat%probe_value(i,mythread) &
                                        + dot_product(d_t, tauab)
          if (dat%do_ground_state) dat%probe_value(i,mythread) = &
            dat%probe_value(i,mythread) + xc%exc(i)*xce%xyzw(i,4)
        end if

!        if (maxval(abs([dsaa,dsbb,dsab,dsba]))<xce%threshold) then
!          d_s = 0
!          f_s = 0
!        end if
!        if (maxval(abs(tauab))<xce%threshold) then
!          d_t = 0
!          f_t = 0
!        end if

        grad_d(i,1,1) = f_r(1)
        grad_p(i,1,1) = d_r(1)

        if (xce%funTyp /= OQP_FUNTYP_LDA) then

          grad_d(i,2:4,1) = &
              (2*f_s(1)+f_s(3)) * drho(1:3,i) &
            + (2*d_s(1)+d_s(3)) * drrho(:,1,i,j)

          grad_p(i,2:4,1) = &
              (2*d_s(1)+d_s(3))*drho(1:3,i)
        end if

        if (xce%funTyp == OQP_FUNTYP_MGGA) then
          grad_d(i,5,1) = f_t(1)
          grad_p(i,5,1) = d_t(1)
        end if

      end do
    end do

    end associate

 end subroutine


!> @brief Compute derivative terms of \sum_ij V^xc_ij P_ij
!>  w.r.t. atomic coordinates
!> @detail this subroutine update terms which should be contracted with
!>  ground state (`grad_d`) and relaxed difference densities (`grad_p`)
!> @note For further details see:
!>  [1] F.Furche, R.Ahlrichs, J.Chem.Phys. 117, p. 7433 (2002)
!>  [2] G.Scalmani et al., J.Chem.Phys. 124, 094107 (2006)
!> @author Vladimir Mironov
 subroutine grad_v_xc(dat, xce, mythread)
    use mod_dft_gridint, only: xc_der1, xc_der2_contr
    class(xc_engine_t) :: xce
    type(xc_consumer_tdg_t) :: dat
    integer :: mythread

    integer :: i, j
    real(kind=fp) :: d_r(2), d_s(3), d_t(2)
    real(kind=fp) :: f_r(2), f_s(3), f_t(2)
    real(kind=fp) :: rhoab(2), tauab(2), sigma(3), dsaa, dsab, dsba, dsbb

    associate ( grad_d  => dat%grad_d(:,:,:,mythread) &
              , grad_p  => dat%grad_p(:,:,:,mythread) &
              , aoG1    => xce%aoG1 &
              , aoV     => xce%aoV  &
              , rrho    => dat%rrho(:,:,:,mythread)  &
              , drrho   => dat%drrho(:,:,:,:,mythread)  &
              , rtau    => dat%rtau(:,:,:,mythread)  &
              , drho    => xce%xclib%drho  &
              , numAOs  => xce%numAOs &
              , numPts  => xce%numPts &
              , nMtx    => dat%nMtx &
              , xc      => xce%XCLib &
              , ids     => xce%XCLib%ids &
              )
    do j = 1, nMtx

      do i = 1, numPts

        rhoab = rrho(1:2,i,j)

        sigma = 0

        tauab = 0

        if (xce%funTyp /= OQP_FUNTYP_LDA) then
          ! Compute sigma terms:
          ! \nabla\rho(D) \dot \nabla\rho(P)
          dsaa = dot_product(drrho(:,1,i,j), drho(1:3,i))
          dsab = dot_product(drrho(:,1,i,j), drho(4:6,i))
          dsbb = dot_product(drrho(:,2,i,j), drho(4:6,i))
          dsba = dot_product(drrho(:,2,i,j), drho(1:3,i))
          sigma = [2*dsaa, 2*dsbb, (dsba+dsab)]
        end if

        if (xce%funTyp == OQP_FUNTYP_MGGA) tauab = rtau(1:2,i,j)

        call xc_der1(xce, xce%hasBeta, i, d_r, d_s, d_t)
        call xc_der2_contr(xce, xce%hasBeta, i, &
                rhoab, sigma, tauab, &
                f_r, f_s, f_t)

        ! Directional XC energy density, already multiplied by the current
        ! quadrature weight because xc_der1 returns the scaled libxc arrays:
        ! q_w = w * delta_P e_xc[D].  It remains separate from the
        ! ground-state XC energy because the moving-grid response is linear
        ! in the relaxed probe P.
        if (dat%do_weight_derivative) then
          if (dat%do_ground_state) then
            dat%probe_value(i,mythread) = &
              xc%exc(i)*xce%xyzw(i,4)
          else
            dat%probe_value(i,mythread) = dot_product(d_r, rhoab)
            if (xce%funTyp /= OQP_FUNTYP_LDA) &
              dat%probe_value(i,mythread) = dat%probe_value(i,mythread) &
                                          + dot_product(d_s, sigma)
            if (xce%funTyp == OQP_FUNTYP_MGGA) &
              dat%probe_value(i,mythread) = dat%probe_value(i,mythread) &
                                          + dot_product(d_t, tauab)
          end if
        end if

!        if (maxval(abs([dsaa,dsbb,dsab,dsba]))<xce%threshold) then
!          d_s = 0
!          f_s = 0
!        end if
!        if (maxval(abs(tauab))<xce%threshold) then
!          d_t = 0
!          f_t = 0
!        end if

        grad_d(i,1,1) = f_r(1)
        grad_p(i,1,1) = d_r(1)

        if (xce%funTyp /= OQP_FUNTYP_LDA) then

          grad_d(i,2:4,1) = &
              2*f_s(1) * drho(1:3,i) &
            +   f_s(3) * drho(4:6,i) &
            + 2*d_s(1) * drrho(:,1,i,j) &
            +   d_s(3) * drrho(:,2,i,j)

          grad_p(i,2:4,1) = &
              2*d_s(1)*drho(1:3,i) &
            +   d_s(3)*drho(4:6,i)
        end if

        if (xce%funTyp == OQP_FUNTYP_MGGA) then
          grad_d(i,5,1) = f_t(1)
          grad_p(i,5,1) = d_t(1)
        end if

        if (.not.xce%hasBeta) cycle

        grad_d(i,1,2) = f_r(2)
        grad_p(i,1,2) = d_r(2)

        if (xce%funTyp /= OQP_FUNTYP_LDA) then
          grad_d(i,2:4,2) = &
              2*f_s(2)*drho(4:6,i) &
            +   f_s(3)*drho(1:3,i) &
            + 2*d_s(2)*drrho(:,2,i,j) &
            +   d_s(3)*drrho(:,1,i,j)

          grad_p(i,2:4,2) = &
              2*d_s(2)*drho(4:6,i) &
            +   d_s(3)*drho(1:3,i)
        end if

        if (xce%funTyp == OQP_FUNTYP_MGGA) then
          grad_d(i,5,2) = f_t(2)
          grad_p(i,5,2) = d_t(2)
        end if

      end do
    end do

    end associate

 end subroutine

!> @brief Compute derivative terms of \sum_ijkl f^xc_ij,kl (X+Y)_ij (X+Y)_kl
!>  w.r.t. atomic coordinates
!> @detail this subroutine update terms which should be contracted with
!>  ground state (`grad_d`) and transition densities (`grad_x`)
!> @note For further details see:
!>  [1] F.Furche, R.Ahlrichs, J.Chem.Phys. 117, p. 7433 (2002)
!>  [2] G.Scalmani et al., J.Chem.Phys. 124, 094107 (2006)
!> @author Vladimir Mironov
 subroutine grad_f_xc_np(dat, xce, mythread)
    use mod_dft_gridint, only: xc_der1, xc_der2_contr, xc_der3_contr
    class(xc_engine_t) :: xce
    type(xc_consumer_tdg_t) :: dat
    integer :: mythread

    integer :: i, j
    real(kind=fp) :: d_r(2), d_s(3), d_t(2)
    real(kind=fp) :: f_r(2), f_s(3), f_t(2)
    real(kind=fp) :: ff_s(3), g_r(2), g_s(3), g_t(2)
    real(kind=fp) :: c(3)
    real(kind=fp) :: rhoab(2), tauab(2), sigma(3), ssigma(3)

    associate ( grad_d  => dat%grad_d(:,:,:,mythread) &
              , grad_x  => dat%grad_x(:,:,:,mythread) &
              , aoG1    => xce%aoG1 &
              , aoV     => xce%aoV  &
              , rrho    => dat%rrho(:,:,:,mythread)  &
              , drrho   => dat%drrho(:,:,:,:,mythread)  &
              , rtau    => dat%rtau(:,:,:,mythread)  &
              , drho    => xce%xclib%drho  &
              , numAOs  => xce%numAOs &
              , numPts  => xce%numPts &
              , nMtx    => dat%nMtx &
              , xc      => xce%XCLib &
              , ids     => xce%XCLib%ids &
              )
    do j = 1, nMtx

      do i = 1, numPts

        rhoab = rrho(1,i,j)

        sigma = 0
        ssigma = 0
        tauab = 0
        if (xce%funTyp /= OQP_FUNTYP_LDA) then
          ! Compute sigma terms for 3rd derivative contraction:
          ! \nabla\rho(D) \dot \nabla\rho(X+Y)
          sigma = 2*dot_product(drrho(:,1,i,j), drho(1:3,i))

          ! \nabla\rho(X+Y) \dot \nabla\rho(X+Y)
          ssigma = 2*dot_product(drrho(:,1,i,j), drrho(:,1,i,j))
        end if

        if (xce%funTyp == OQP_FUNTYP_MGGA) tauab = rtau(1,i,j)

        call xc_der1(xce, .false., i, d_r, d_s, d_t)
        call xc_der2_contr(xce, .false., i, &
                rhoab, sigma, tauab, &
                f_r, f_s, f_t)

        call xc_der3_contr(xce, i, &
                rhoab, sigma, tauab, &
                ssigma, &
                ff_s, &
                g_r, g_s, g_t)

        if (dat%do_weight_derivative) then
          dat%probe_value(i,mythread) = dat%probe_value(i,mythread) &
                                      + dot_product(f_r,rhoab)
          if (xce%funTyp /= OQP_FUNTYP_LDA) &
            dat%probe_value(i,mythread) = dat%probe_value(i,mythread) &
                                        + dot_product(f_s,sigma)
          if (xce%funTyp == OQP_FUNTYP_MGGA) &
            dat%probe_value(i,mythread) = dat%probe_value(i,mythread) &
                                        + dot_product(f_t,tauab)
        end if

!        if (maxval(abs([dsaa,dsbb,dsab,dsba]))<xce%threshold) then
!          f_s = 0
!          g_s = 0
!          ff_s = 0
!        end if
!        if (maxval(abs(tauab))<xce%threshold) then
!          f_t = 0
!          g_t = 0
!        end if

        grad_x(i,1,1) = 2*f_r(1)
        grad_d(i,1,1) = grad_d(i,1,1) + g_r(1)

        if (xce%funTyp /= OQP_FUNTYP_LDA) then

          c = &
              (2*f_s(1)+f_s(3))*drho(1:3,i) &
            + (2*d_s(1)+d_s(3))*drrho(:,1,i,j)

          grad_x(i,2:4,1) = 2*c

          grad_d(i,2:4,1) = grad_d(i,2:4,1) &
            + 2*(2*f_s(1)+f_s(3))*drrho(:,1,i,j) &
            +   (2*g_s(1)+g_s(3)) * drho(1:3,i)

        end if

        if (xce%funTyp == OQP_FUNTYP_MGGA) then
           grad_x(i,5,1) = 2*f_t(1)
           grad_d(i,5,1) = grad_d(i,5,1) + g_t(1)
        end if

      end do
    end do

    end associate

 end subroutine

!> @brief Compute derivative terms of \sum_ijkl f^xc_ij,kl (X+Y)_ij (X+Y)_kl
!>  w.r.t. atomic coordinates
!> @detail this subroutine update terms which should be contracted with
!>  ground state (`grad_d`) and transition densities (`grad_x`)
!> @note For further details see:
!>  [1] F.Furche, R.Ahlrichs, J.Chem.Phys. 117, p. 7433 (2002)
!>  [2] G.Scalmani et al., J.Chem.Phys. 124, 094107 (2006)
!> @author Vladimir Mironov
 subroutine grad_f_xc(dat, xce, mythread)
    use mod_dft_gridint, only: xc_der1, xc_der2_contr, xc_der3_contr
    class(xc_engine_t) :: xce
    type(xc_consumer_tdg_t) :: dat
    integer :: mythread

    integer :: i, j
    real(kind=fp) :: d_r(2), d_s(3), d_t(2)
    real(kind=fp) :: f_r(2), f_s(3), f_t(2)
    real(kind=fp) :: ff_s(3), g_r(2), g_s(3), g_t(2)
    real(kind=fp) :: c(3)
    real(kind=fp) :: rhoab(2), tauab(2), sigma(3), ssigma(3), dsaa, dsab, dsba, dsbb
    real(kind=fp) :: ssaa, ssab, ssbb

    associate ( grad_d  => dat%grad_d(:,:,:,mythread) &
              , grad_x  => dat%grad_x(:,:,:,mythread) &
              , aoG1    => xce%aoG1 &
              , aoV     => xce%aoV  &
              , rrho    => dat%rrho(:,:,:,mythread)  &
              , drrho   => dat%drrho(:,:,:,:,mythread)  &
              , rtau    => dat%rtau(:,:,:,mythread)  &
              , drho    => xce%xclib%drho  &
              , numAOs  => xce%numAOs &
              , numPts  => xce%numPts &
              , nMtx    => dat%nMtx &
              , xc      => xce%XCLib &
              , ids     => xce%XCLib%ids &
              )
    do j = 1, nMtx

      do i = 1, numPts

        rhoab = rrho(1:2,i,j)

        sigma = 0
        ssigma = 0

        tauab = 0

        if (xce%funTyp /= OQP_FUNTYP_LDA) then
          ! Compute sigma terms for 3rd derivative contraction:
          ! \nabla\rho(D) \dot \nabla\rho(X+Y)
          dsaa = dot_product(drrho(:,1,i,j), drho(1:3,i))
          dsab = dot_product(drrho(:,1,i,j), drho(4:6,i))
          dsbb = dot_product(drrho(:,2,i,j), drho(4:6,i))
          dsba = dot_product(drrho(:,2,i,j), drho(1:3,i))

          ! \nabla\rho(X+Y) \dot \nabla\rho(X+Y)
          ssaa = dot_product(drrho(:,1,i,j), drrho(:,1,i,j))
          ssab = dot_product(drrho(:,1,i,j), drrho(:,2,i,j))
          ssbb = dot_product(drrho(:,2,i,j), drrho(:,2,i,j))

          ! Compute \sigma_aa, \sigma_bb, \sigma_ab for the above:
          sigma = [2*dsaa, 2*dsbb, (dsab+dsba)]
          ssigma = [2*ssaa, 2*ssbb, 2*ssab]
        end if

        if (xce%funTyp == OQP_FUNTYP_MGGA) tauab = rtau(1:2,i,j)

        call xc_der1(xce, xce%hasBeta, i, d_r, d_s, d_t)
        call xc_der2_contr(xce, xce%hasBeta, i, &
                rhoab, sigma, tauab, &
                f_r, f_s, f_t)

        call xc_der3_contr(xce, i, &
                rhoab, sigma, tauab, &
                ssigma, &
                ff_s, &
                g_r, g_s, g_t)

!        if (maxval(abs([dsaa,dsbb,dsab,dsba]))<xce%threshold) then
!          f_s = 0
!          g_s = 0
!          ff_s = 0
!        end if
!        if (maxval(abs(tauab))<xce%threshold) then
!          f_t = 0
!          g_t = 0
!        end if

        grad_x(i,1,1) = 2*f_r(1)
        grad_d(i,1,1) = grad_d(i,1,1) + g_r(1)

        if (xce%funTyp /= OQP_FUNTYP_LDA) then

          c = &
              2*f_s(1)*drho(1:3,i) &
            +   f_s(3)*drho(4:6,i) &
            + 2*d_s(1)*drrho(:,1,i,j) &
            +   d_s(3)*drrho(:,2,i,j)

          grad_x(i,2:4,1) = 2*c

          c = &
            + 2*f_s(1) * drrho(:,1,i,j) &
            +   f_s(3) * drrho(:,2,i,j)
          grad_d(i,2:4,1) = grad_d(i,2:4,1) &
            + 2*c &
            + 2*g_s(1) * drho(1:3,i) &
            +   g_s(3) * drho(4:6,i)

        end if

        if (xce%funTyp == OQP_FUNTYP_MGGA) then
           grad_x(i,5,1) = 2*f_t(1)
           grad_d(i,5,1) = grad_d(i,5,1) + g_t(1)
        end if

        if (.not.xce%hasBeta) cycle

        grad_x(i,1,2) = 2*f_r(2)
        grad_d(i,1,2) = grad_d(i,1,2) + g_r(2)

        if (xce%funTyp /= OQP_FUNTYP_LDA) then

          c = &
              2*f_s(2)*drho(4:6,i) &
            +   f_s(3)*drho(1:3,i) &
            + 2*d_s(2)*drrho(:,2,i,j) &
            +   d_s(3)*drrho(:,1,i,j)

          grad_x(i,2:4,2) = 2*c

          c = &
            + 2*f_s(2) * drrho(:,2,i,j) &
            +   f_s(3) * drrho(:,1,i,j)

          grad_d(i,2:4,2) = grad_d(i,2:4,2) &
            + 2*c &
            + 2*g_s(2) * drho(4:6,i) &
            +   g_s(3) * drho(1:3,i)

        end if

        if (xce%funTyp == OQP_FUNTYP_MGGA) then
           grad_x(i,5,2) = 2*f_t(2)
           grad_d(i,5,2) = grad_d(i,5,2) + g_t(2)
        end if

      end do
    end do

    end associate

 end subroutine

!> @brief Compute derivative XC contribution to the TD-DFT KS-like matrices
!> @param[in]    basis     basis set
!> @param[in]    wf        density matrix/orbitals
!> @param[inout] fx        fock-like matrices
!> @param[inout] dx        densities
!> @param[in]    nMtx      number of density/Fock-like matrices
!> @param[in]    threshold tolerance
!> @param[in]    isGGA     .TRUE. if GGA/mGGA functional used
!> @param[in]    infos     OQP metadata
!> @author Vladimir Mironov
  subroutine utddft_xc_gradient(basis, molGrid, dedft, &
                  da, db, pa, pb, xa, xb, &
                  nMtx, threshold, infos, &
                  include_ground_state, include_weight_derivative, &
                  weight_derivative_only)
!$  use omp_lib, only: omp_get_num_threads, omp_get_thread_num
    use basis_tools, only: basis_set
    use mod_dft_gridint, only: xc_options_t, run_xc
    use types, only: information
    use mod_dft_molgrid, only: dft_grid_t
    use messages, only: show_message, with_abort

    implicit none

    type(information), target, intent(in) :: infos
    type(dft_grid_t), target, intent(in) :: molGrid
    real(kind=fp), intent(inout) :: dedft(:,:)

    type(basis_set) :: basis
    integer, intent(in) :: nMtx
    real(kind=fp), intent(inout), contiguous, target :: da(:,:), db(:,:)
    real(kind=fp), intent(inout), target :: pa(:,:,:), pb(:,:,:)
    real(kind=fp), intent(inout), optional, target :: &
            xa(:,:,:), xb(:,:,:)
    real(kind=fp), intent(in) :: threshold
    logical, intent(in), optional :: include_ground_state
    logical, intent(in), optional :: include_weight_derivative
    logical, intent(in), optional :: weight_derivative_only

    type(xc_consumer_tdg_t) :: dat
    type(xc_options_t) :: xc_opts

    integer :: i, j, nbf, nxcder
    logical :: doFxc, requested_weight_derivative, requested_weight_only

    nbf = ubound(da,1)
    doFxc = present(xa)
    requested_weight_derivative = .false.
    if (present(include_weight_derivative)) &
      requested_weight_derivative = include_weight_derivative
    requested_weight_only = .false.
    if (present(weight_derivative_only)) &
      requested_weight_only = weight_derivative_only
    if (requested_weight_derivative .and. doFxc) &
      call show_message('XC moving-grid response is available only for '// &
                        'the linear-probe branch (xa/xb absent).', with_abort)
    if (requested_weight_derivative .and. nMtx /= 1) &
      call show_message('XC moving-grid response requires nMtx=1.', with_abort)
    if (requested_weight_only .and. .not. requested_weight_derivative) &
      call show_message('XC weight_derivative_only requires '// &
                        'include_weight_derivative=.true..', with_abort)

    ! Scale densities by B.F. norms
      do i = 1, nbf
        da(:,i) = da(:,i) &
                 * basis%bfnrm(i) &
                 * basis%bfnrm(:)
        db(:,i) = db(:,i) &
                 * basis%bfnrm(i) &
                 * basis%bfnrm(:)
      end do
    do j = 1, nMtx
      do i = 1, nbf
        pa(:,i,j) = pa(:,i,j) &
                 * basis%bfnrm(i) &
                 * basis%bfnrm(:)
        pb(:,i,j) = pb(:,i,j) &
                 * basis%bfnrm(i) &
                 * basis%bfnrm(:)
      end do
    end do

    if (doFxc) then
      do j = 1, nMtx
        do i = 1, nbf
          xa(:,i,j) = xa(:,i,j) &
                   * basis%bfnrm(i) &
                   * basis%bfnrm(:)
          xb(:,i,j) = xb(:,i,j) &
                   * basis%bfnrm(i) &
                   * basis%bfnrm(:)
        end do
      end do
    end if

    nxcder = 2
    if (doFxc) nxcder = 3

    xc_opts%isGGA = infos%functional%needGrd
    xc_opts%needTau = infos%functional%needTau
    xc_opts%hasBeta = .true.
    xc_opts%isWFVecs = .false.
    xc_opts%numAOs = nbf
    xc_opts%maxPts = molGrid%maxSlicePts
    xc_opts%limPts = molGrid%maxNRadTimesNAng
    xc_opts%numAtoms = infos%mol_prop%natom
    xc_opts%maxAngMom = basis%mxam
    xc_opts%functional => infos%functional
    xc_opts%nDer = 1
    xc_opts%nXCDer = nxcder
    xc_opts%numOccAlpha = infos%mol_prop%nelec_A
    xc_opts%numOccBeta = infos%mol_prop%nelec_B
    xc_opts%wfAlpha => da
    xc_opts%wfBeta  => db
    xc_opts%dft_threshold = threshold
    xc_opts%molGrid => molGrid

    dat%pa => pa
    dat%pb => pb
    if (doFxc) then
      dat%xa => xa
      dat%xb => xb
    end if
    dat%nMtx = nMtx
    dat%do_fxc = doFxc
    dat%do_ground_state = .true.
    if (present(include_ground_state)) &
      dat%do_ground_state = include_ground_state
    dat%do_weight_derivative = requested_weight_derivative
    if (dat%do_weight_derivative) then
      dat%part_fun_type = molGrid%partFunType
      dat%has_surface_shift = molGrid%hasSurfaceShift
      dat%atom_xyz = infos%atoms%xyz
      dat%dummy_atom = molGrid%dummyAtom
      dat%surface_shift = molGrid%surfaceShift
    end if

    call dat%pe%init(infos%mpiinfo%comm, infos%mpiinfo%usempi)

    call run_xc(xc_opts, dat, basis)

    ! Scale densities back
    do i = 1, nbf
      da(:,i) = da(:,i) &
             / basis%bfnrm(i) &
             / basis%bfnrm(:)
      db(:,i) = db(:,i) &
               / basis%bfnrm(i) &
               / basis%bfnrm(:)
    end do
    do j = 1, nMtx
      do i = 1, nbf
        pa(:,i,j) = pa(:,i,j) &
                 / basis%bfnrm(i) &
                 / basis%bfnrm(:)
        pb(:,i,j) = pb(:,i,j) &
                 / basis%bfnrm(i) &
                 / basis%bfnrm(:)
      end do
    end do

    if (doFxc) then
      do j = 1, nMtx
        do i = 1, nbf
          xa(:,i,j) = xa(:,i,j) &
                   / basis%bfnrm(i) &
                   / basis%bfnrm(:)
          xb(:,i,j) = xb(:,i,j) &
                   / basis%bfnrm(i) &
                   / basis%bfnrm(:)
        end do
      end do
    end if

    if (.not. requested_weight_only) then
      do j = 1, basis%nshell
        associate (atom => basis%origin(j), &
                   offset => basis%ao_offset(j), &
                   naos => basis%naos(j))
          dedft(1, atom) = dedft(1, atom)-sum(dat%bfGrad(offset:offset+naos-1, 1, 1))
          dedft(2, atom) = dedft(2, atom)-sum(dat%bfGrad(offset:offset+naos-1, 2, 1))
          dedft(3, atom) = dedft(3, atom)-sum(dat%bfGrad(offset:offset+naos-1, 3, 1))
        end associate
      end do
    end if

    if (dat%do_weight_derivative) dedft = dedft + dat%nucGrad(:,:,1)

    call dat%clean()
  end subroutine

!-------------------------------------------------------------------------------

!> @brief Compute derivative XC contribution to the TD-DFT KS-like matrices
!> @param[in]    basis     basis set
!> @param[in]    wf        density matrix/orbitals
!> @param[inout] fx        fock-like matrices
!> @param[inout] dx        densities
!> @param[in]    nMtx      number of density/Fock-like matrices
!> @param[in]    threshold tolerance
!> @param[in]    isGGA     .TRUE. if GGA/mGGA functional used
!> @param[in]    infos     OQP metadata
!> @author Vladimir Mironov
  subroutine tddft_xc_gradient(basis, molGrid, dedft, &
                  da, pa, xa, &
                  nMtx, threshold, infos, include_weight_derivative)
!$  use omp_lib, only: omp_get_num_threads, omp_get_thread_num
    use basis_tools, only: basis_set
    use mod_dft_gridint, only: xc_options_t, run_xc
    use types, only: information
    use mod_dft_molgrid, only: dft_grid_t

    implicit none

    type(information), target, intent(in) :: infos
    type(dft_grid_t), target, intent(in) :: molGrid
    real(kind=fp), intent(inout) :: dedft(:,:)

    type(basis_set) :: basis
    integer, intent(in) :: nMtx
    real(kind=fp), intent(inout), contiguous, target :: da(:,:)
    real(kind=fp), intent(inout), target :: pa(:,:,:)
    real(kind=fp), intent(inout), optional, target :: xa(:,:,:)
    real(kind=fp), intent(in) :: threshold
    logical, intent(in), optional :: include_weight_derivative

    type(xc_consumer_tdg_t) :: dat
    type(xc_options_t) :: xc_opts

    integer :: i, j, nbf, nxcder
    logical :: doFxc, doWeight

    nbf = ubound(da,1)
    ! Scale densities by B.F. norms
      do i = 1, nbf
        da(:,i) = da(:,i) &
                 * basis%bfnrm(i) &
                 * basis%bfnrm(:)
      end do
    do j = 1, nMtx
      do i = 1, nbf
        pa(:,i,j) = pa(:,i,j) &
                 * basis%bfnrm(i) &
                 * basis%bfnrm(:)
      end do
    end do

    doFxc = present(xa)
    doWeight = .false.
    if (present(include_weight_derivative)) doWeight = include_weight_derivative

    if (doFxc) then
      do j = 1, nMtx
        do i = 1, nbf
          xa(:,i,j) = xa(:,i,j) &
                   * basis%bfnrm(i) &
                   * basis%bfnrm(:)
        end do
      end do
    end if

    nxcder = 2
    if (doFxc) nxcder = 3

    xc_opts%isGGA = infos%functional%needGrd
    xc_opts%needTau = infos%functional%needTau
    xc_opts%hasBeta = .false.
    xc_opts%isWFVecs = .false.
    xc_opts%numAOs = nbf
    xc_opts%maxPts = molGrid%maxSlicePts
    xc_opts%limPts = molGrid%maxNRadTimesNAng
    xc_opts%numAtoms = infos%mol_prop%natom
    xc_opts%maxAngMom = basis%mxam
    xc_opts%functional => infos%functional
    xc_opts%nDer = 1
    xc_opts%nXCDer = nxcder
    xc_opts%numOccAlpha = infos%mol_prop%nelec_A
    xc_opts%numOccBeta = infos%mol_prop%nelec_B
    xc_opts%wfAlpha => da
    xc_opts%dft_threshold = threshold
    xc_opts%molGrid => molGrid

    dat%pa => pa
    if (doFxc) then
      dat%xa => xa
    end if
    dat%nMtx = nMtx
    dat%do_fxc = doFxc
    dat%do_weight_derivative = doWeight
    if (doWeight) then
      dat%part_fun_type = molGrid%partFunType
      dat%has_surface_shift = molGrid%hasSurfaceShift
      dat%atom_xyz = infos%atoms%xyz
      dat%dummy_atom = molGrid%dummyAtom
      dat%surface_shift = molGrid%surfaceShift
    end if

    call dat%pe%init(infos%mpiinfo%comm, infos%mpiinfo%usempi)

    call run_xc(xc_opts, dat, basis)

    ! Scale densities back
    do i = 1, nbf
      da(:,i) = da(:,i) &
             / basis%bfnrm(i) &
             / basis%bfnrm(:)
    end do
    do j = 1, nMtx
      do i = 1, nbf
        pa(:,i,j) = pa(:,i,j) &
                 / basis%bfnrm(i) &
                 / basis%bfnrm(:)
      end do
    end do

    if (doFxc) then
      do j = 1, nMtx
        do i = 1, nbf
          xa(:,i,j) = xa(:,i,j) &
                   / basis%bfnrm(i) &
                   / basis%bfnrm(:)
        end do
      end do
    end if

    ! Factor 2 is because only alpha contribution to gradient is computed abouve
    ! beta contribution is equal to alpha in RHF case
    do j = 1, basis%nshell
      associate (atom => basis%origin(j), &
                 offset => basis%ao_offset(j), &
                 naos => basis%naos(j))
        dedft(1, atom) = dedft(1, atom)-2*sum(dat%bfGrad(offset:offset+naos-1, 1, 1))
        dedft(2, atom) = dedft(2, atom)-2*sum(dat%bfGrad(offset:offset+naos-1, 2, 1))
        dedft(3, atom) = dedft(3, atom)-2*sum(dat%bfGrad(offset:offset+naos-1, 3, 1))
      end associate
    end do

    if (dat%do_weight_derivative) dedft = dedft + dat%nucGrad(:,:,1)

    call dat%clean()
  end subroutine

end module mod_dft_gridint_tdxc_grad
