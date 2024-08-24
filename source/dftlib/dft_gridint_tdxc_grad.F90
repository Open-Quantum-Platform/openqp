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
    real(kind=fp), pointer :: pa(:,:,:)
    real(kind=fp), pointer :: pb(:,:,:)
    real(kind=fp), pointer :: xa(:,:,:)
    real(kind=fp), pointer :: xb(:,:,:)
    real(kind=fp), allocatable :: rrho(:,:,:,:)
    real(kind=fp), allocatable :: drrho(:,:,:,:,:)
    real(kind=fp), allocatable :: rtau(:,:,:,:)
    real(kind=fp), allocatable :: bfgrad(:,:,:)
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
    integer :: nspin, nterms, nDeriv
    call self%clean()
    nterms = 1
    if (xce%funTyp /= OQP_FUNTYP_LDA) nterms = nterms + 3
    if (xce%funTyp == OQP_FUNTYP_MGGA) nterms = nterms + 1

    nspin = merge(2, 1, xce%hasBeta)
    nDeriv = merge(2, 1, self%do_fxc)
    allocate( &
        self%bfgrad(xce%numAOs, 3, nthreads) &
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
  end subroutine

!-------------------------------------------------------------------------------

  subroutine parallel_stop(self)
    implicit none
    class(xc_consumer_tdg_t), intent(inout) :: self
    if (ubound(self%bfGrad,3) <= 1) return
    self%bfGrad(:,:,lbound(self%bfGrad,3)) = sum(self%bfGrad, dim=3)
  end subroutine

!-------------------------------------------------------------------------------

  subroutine clean(self)
    implicit none
    class(xc_consumer_tdg_t), intent(inout) :: self
    if (allocated(self%bfgrad)) deallocate(self%bfgrad)
    if (allocated(self%rrho)) deallocate(self%rrho)
    if (allocated(self%drrho)) deallocate(self%drrho)
    if (allocated(self%rtau)) deallocate(self%rtau)
    if (allocated(self%tmpGrad_)) deallocate(self%tmpGrad_)
    if (allocated(self%tmp_)) deallocate(self%tmp_)
    if (allocated(self%tmpV_)) deallocate(self%tmpV_)
    if (allocated(self%tmpG1_)) deallocate(self%tmpG1_)
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

   end associate

 end subroutine

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
                  nMtx, threshold, infos)
!$  use omp_lib, only: omp_get_num_threads, omp_get_thread_num
    use basis_tools, only: basis_set
    use mod_dft_gridint, only: xc_options_t, run_xc
    use types, only: information
    use mod_dft_molgrid, only: dft_grid_t

    implicit none

    type(information), target, intent(in) :: infos
    type(dft_grid_t), target, intent(in) :: molGrid
    real(kind=fp), intent(out) :: dedft(:,:)

    type(basis_set) :: basis
    integer, intent(in) :: nMtx
    real(kind=fp), intent(inout), contiguous, target :: da(:,:), db(:,:)
    real(kind=fp), intent(inout), target :: pa(:,:,:), pb(:,:,:)
    real(kind=fp), intent(inout), optional, target :: &
            xa(:,:,:), xb(:,:,:)
    real(kind=fp), intent(in) :: threshold

    type(xc_consumer_tdg_t) :: dat
    type(xc_options_t) :: xc_opts

    integer :: i, j, nbf, nxcder
    logical :: doFxc

    nbf = ubound(da,1)

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

    doFxc = present(xa)

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
    xc_opts%maxAngMom = basis%basis_max_angular_momentum
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

    do j = 1, basis%nshell
      associate (kat => basis%katom(j), &
                 kloc => basis%kloc(j), &
                 knbf => basis%kmax(j)-basis%kmin(j)+1)
        dedft(1, kat) = dedft(1, kat)-sum(dat%bfGrad(kloc:kloc+knbf-1, 1, 1))
        dedft(2, kat) = dedft(2, kat)-sum(dat%bfGrad(kloc:kloc+knbf-1, 2, 1))
        dedft(3, kat) = dedft(3, kat)-sum(dat%bfGrad(kloc:kloc+knbf-1, 3, 1))
      end associate
    end do

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
                  nMtx, threshold, infos)
!$  use omp_lib, only: omp_get_num_threads, omp_get_thread_num
    use basis_tools, only: basis_set
    use mod_dft_gridint, only: xc_options_t, run_xc
    use types, only: information
    use mod_dft_molgrid, only: dft_grid_t

    implicit none

    type(information), target, intent(in) :: infos
    type(dft_grid_t), target, intent(in) :: molGrid
    real(kind=fp), intent(out) :: dedft(:,:)

    type(basis_set) :: basis
    integer, intent(in) :: nMtx
    real(kind=fp), intent(inout), contiguous, target :: da(:,:)
    real(kind=fp), intent(inout), target :: pa(:,:,:)
    real(kind=fp), intent(inout), optional, target :: xa(:,:,:)
    real(kind=fp), intent(in) :: threshold

    type(xc_consumer_tdg_t) :: dat
    type(xc_options_t) :: xc_opts

    integer :: i, j, nbf, nxcder
    logical :: doFxc

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
    xc_opts%maxAngMom = basis%basis_max_angular_momentum
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
      associate (kat => basis%katom(j), &
                 kloc => basis%kloc(j), &
                 knbf => basis%kmax(j)-basis%kmin(j)+1)
        dedft(1, kat) = dedft(1, kat)-2*sum(dat%bfGrad(kloc:kloc+knbf-1, 1, 1))
        dedft(2, kat) = dedft(2, kat)-2*sum(dat%bfGrad(kloc:kloc+knbf-1, 2, 1))
        dedft(3, kat) = dedft(3, kat)-2*sum(dat%bfGrad(kloc:kloc+knbf-1, 3, 1))
      end associate
    end do

    call dat%clean()
  end subroutine

end module mod_dft_gridint_tdxc_grad
