module mod_dft_gridint_giao

  use precision, only: fp
  use mod_dft_gridint, only: xc_engine_t, xc_consumer_t, OQP_FUNTYP_LDA
  use oqp_linalg

  implicit none

!-------------------------------------------------------------------------------
! London (GIAO) derivative of the exchange-correlation potential (the
! "vxc_giao" term) in the first-order magnetic Hamiltonian for GIAO NMR.
! Per spin, accumulates three real matrices V(:,:,t); the caller antisymmetrizes
! (V - V^T) and subtracts from h1.
!   V_t[mu,nu] += aow_mu * ig_t,nu  +  aoV_mu * sum_g wv_g ipig_{g,t},nu
!   ig_t,nu       = 0.5 (R_nu x r)_t aoV_nu                       (London AO value derivative)
!   ipig_{g,t},nu = 0.5[(R_nu x r)_t aoG1_nu,g + (R_nu x e_g)_t aoV_nu] (London AO gradient derivative)
!   aow_mu        = vrho*aoV_mu + sum_g wv_g aoG1_mu,g
! NOTE: d1dr/d1ds already carry the grid weight (XCLib%compute(.,wts)).
! NOTE: when xce%skip_p (no AO pruning) numAOs_p==numAOs and the AOs are in
!   natural order, so the full AO index is k itself; indices_p(k) is only valid
!   for the pruned (skip_p == .false.) case.
!-------------------------------------------------------------------------------

  type, extends(xc_consumer_t) :: xc_consumer_giao_t
    real(kind=fp), allocatable :: aocent(:,:)
    real(kind=fp), allocatable :: va2(:,:), vb2(:,:)
    real(kind=fp), allocatable :: vmat_(:,:), igs_(:,:), aow_(:,:)
  contains
    procedure :: parallel_start
    procedure :: parallel_stop
    procedure :: update
    procedure :: postUpdate
    procedure :: clean
    procedure :: getPtr
  end type

  private
  public xc_consumer_giao_t
  public giao_vxc

contains

  subroutine parallel_start(self, xce, nthreads)
    class(xc_consumer_giao_t), target, intent(inout) :: self
    class(xc_engine_t), intent(in) :: xce
    integer, intent(in) :: nthreads
    integer :: nSpin
    nSpin = 1
    if (xce%hasBeta) nSpin = 2
    allocate(self%va2(xce%numAOs*xce%numAOs*3, nthreads), source=0.0_fp)
    if (xce%hasBeta) allocate(self%vb2(xce%numAOs*xce%numAOs*3, nthreads), source=0.0_fp)
    allocate(self%vmat_(xce%numAOs*xce%numAOs*3*nSpin, nthreads), source=0.0_fp)
    allocate(self%igs_(xce%numAOs*xce%maxPts*3, nthreads), source=0.0_fp)
    allocate(self%aow_(xce%numAOs*xce%maxPts, nthreads), source=0.0_fp)
  end subroutine

  subroutine parallel_stop(self)
    class(xc_consumer_giao_t), intent(inout) :: self
    if (ubound(self%va2,2) /= 1) self%va2(:,1) = sum(self%va2, dim=2)
    call self%pe%allreduce(self%va2(:,1), size(self%va2(:,1)))
    if (allocated(self%vb2)) then
      if (ubound(self%vb2,2) /= 1) self%vb2(:,1) = sum(self%vb2, dim=2)
      call self%pe%allreduce(self%vb2(:,1), size(self%vb2(:,1)))
    end if
  end subroutine

  subroutine clean(self)
    class(xc_consumer_giao_t), intent(inout) :: self
    if (allocated(self%va2)) deallocate(self%va2)
    if (allocated(self%vb2)) deallocate(self%vb2)
    if (allocated(self%vmat_)) deallocate(self%vmat_)
    if (allocated(self%igs_)) deallocate(self%igs_)
    if (allocated(self%aow_)) deallocate(self%aow_)
    if (allocated(self%aocent)) deallocate(self%aocent)
  end subroutine

  subroutine getPtr(self, nAOp, nPts, nSpin, mythread, vmat, ig, aow)
    class(xc_consumer_giao_t), target, intent(inout) :: self
    integer, intent(in) :: nAOp, nPts, nSpin, mythread
    real(kind=fp), pointer, intent(out) :: vmat(:,:,:,:), ig(:,:,:), aow(:,:)
    vmat(1:nAOp,1:nAOp,1:3,1:nSpin) => self%vmat_(1:nAOp*nAOp*3*nSpin, mythread)
    ig(1:nAOp,1:nPts,1:3) => self%igs_(1:nAOp*nPts*3, mythread)
    aow(1:nAOp,1:nPts) => self%aow_(1:nAOp*nPts, mythread)
  end subroutine

  subroutine update(self, xce, mythread)
    class(xc_consumer_giao_t), intent(inout) :: self
    class(xc_engine_t), intent(in) :: xce
    integer :: mythread
    real(kind=fp), pointer :: vmat(:,:,:,:), ig(:,:,:), aow(:,:)
    real(kind=fp), allocatable :: g2(:,:), rc(:,:)
    integer :: nSpin, sp, i, k, t
    real(kind=fp) :: R(3), r3(3), wv(3), cr(3), cre(3,3)

    nSpin = 1
    if (xce%hasBeta) nSpin = 2

    associate( nAOp => xce%numAOs_p, nPts => xce%numPts, &
               aoV => xce%aoV, aoG1 => xce%aoG1, &
               d1dr => xce%XCLib%d1dr, d1ds => xce%XCLib%d1ds, &
               drho => xce%XCLib%drho, ra => xce%XCLib%ids%ra, &
               rb => xce%XCLib%ids%rb, ga => xce%XCLib%ids%ga, &
               gb => xce%XCLib%ids%gb, gc => xce%XCLib%ids%gc, &
               indices => xce%indices_p )

      call self%getPtr(nAOp, nPts, nSpin, mythread, vmat, ig, aow)
      vmat = 0.0_fp

      ! Cache the AO centers in pruned order (skip_p -> natural order k).
      allocate(rc(3,nAOp))
      do k = 1, nAOp
        if (xce%skip_p) then
          rc(:,k) = self%aocent(:, k)
        else
          rc(:,k) = self%aocent(:, indices(k))
        end if
      end do

      ! GIAO value weight ig_t,nu = 0.5 (R_nu x r)_t aoV_nu
      do i = 1, nPts
        r3 = xce%xyzw(i,1:3)
        do k = 1, nAOp
          R = rc(:,k)
          cr(1) = R(2)*r3(3) - R(3)*r3(2)
          cr(2) = R(3)*r3(1) - R(1)*r3(3)
          cr(3) = R(1)*r3(2) - R(2)*r3(1)
          ig(k,i,1) = 0.5_fp*cr(1)*aoV(k,i)
          ig(k,i,2) = 0.5_fp*cr(2)*aoV(k,i)
          ig(k,i,3) = 0.5_fp*cr(3)*aoV(k,i)
        end do
      end do

      if (xce%funTyp /= OQP_FUNTYP_LDA) allocate(g2(nAOp,nPts))

      do sp = 1, nSpin
        do i = 1, nPts
          if (sp == 1) then
            aow(:,i) = d1dr(ra,i)*aoV(:,i)
          else
            aow(:,i) = d1dr(rb,i)*aoV(:,i)
          end if
        end do
        if (xce%funTyp /= OQP_FUNTYP_LDA) then
          do i = 1, nPts
            if (sp == 1) then
              wv = 2.0_fp*d1ds(ga,i)*drho(1:3,i) + d1ds(gc,i)*drho(4:6,i)
            else
              wv = 2.0_fp*d1ds(gb,i)*drho(4:6,i) + d1ds(gc,i)*drho(1:3,i)
            end if
            aow(:,i) = aow(:,i) + wv(1)*aoG1(:,i,1) + wv(2)*aoG1(:,i,2) + wv(3)*aoG1(:,i,3)
          end do
        end if

        do t = 1, 3
          call dgemm('N','T', nAOp, nAOp, nPts, 1.0_fp, &
                     aow, nAOp, ig(:,:,t), nAOp, 1.0_fp, vmat(:,:,t,sp), nAOp)
        end do

        if (xce%funTyp /= OQP_FUNTYP_LDA) then
          do t = 1, 3
            do i = 1, nPts
              r3 = xce%xyzw(i,1:3)
              if (sp == 1) then
                wv = 2.0_fp*d1ds(ga,i)*drho(1:3,i) + d1ds(gc,i)*drho(4:6,i)
              else
                wv = 2.0_fp*d1ds(gb,i)*drho(4:6,i) + d1ds(gc,i)*drho(1:3,i)
              end if
              do k = 1, nAOp
                R = rc(:,k)
                cr(1) = R(2)*r3(3) - R(3)*r3(2)
                cr(2) = R(3)*r3(1) - R(1)*r3(3)
                cr(3) = R(1)*r3(2) - R(2)*r3(1)
                cre(1,1)=0.0_fp; cre(2,1)= R(3);  cre(3,1)=-R(2)
                cre(1,2)=-R(3);  cre(2,2)=0.0_fp; cre(3,2)= R(1)
                cre(1,3)= R(2);  cre(2,3)=-R(1);  cre(3,3)=0.0_fp
                g2(k,i) = 0.5_fp*( wv(1)*(cr(t)*aoG1(k,i,1)+cre(t,1)*aoV(k,i)) &
                                 + wv(2)*(cr(t)*aoG1(k,i,2)+cre(t,2)*aoV(k,i)) &
                                 + wv(3)*(cr(t)*aoG1(k,i,3)+cre(t,3)*aoV(k,i)) )
              end do
            end do
            call dgemm('N','T', nAOp, nAOp, nPts, 1.0_fp, &
                       aoV, nAOp, g2, nAOp, 1.0_fp, vmat(:,:,t,sp), nAOp)
          end do
        end if
      end do
      if (allocated(g2)) deallocate(g2)
      deallocate(rc)
    end associate
  end subroutine

  subroutine postUpdate(self, xce, mythread)
    class(xc_consumer_giao_t), intent(inout) :: self
    class(xc_engine_t), intent(in) :: xce
    integer :: mythread
    real(kind=fp), pointer :: vmat(:,:,:,:), ig(:,:,:), aow(:,:)
    real(kind=fp), pointer :: va(:,:,:), vb(:,:,:)
    integer :: nSpin, t, nAO

    nSpin = 1
    if (xce%hasBeta) nSpin = 2
    nAO = xce%numAOs
    associate( nAOp => xce%numAOs_p, indices => xce%indices_p )
      call self%getPtr(nAOp, max(xce%numPts,1), nSpin, mythread, vmat, ig, aow)
      call mapfull(self%va2(:,mythread), nAO, va)
      if (xce%hasBeta) call mapfull(self%vb2(:,mythread), nAO, vb)
      if (xce%skip_p) then
        do t = 1, 3
          va(:,:,t) = va(:,:,t) + vmat(:,:,t,1)
          if (xce%hasBeta) vb(:,:,t) = vb(:,:,t) + vmat(:,:,t,2)
        end do
      else
        do t = 1, 3
          va(indices(1:nAOp), indices(1:nAOp), t) = &
            va(indices(1:nAOp), indices(1:nAOp), t) + vmat(:,:,t,1)
          if (xce%hasBeta) &
            vb(indices(1:nAOp), indices(1:nAOp), t) = &
              vb(indices(1:nAOp), indices(1:nAOp), t) + vmat(:,:,t,2)
        end do
      end if
    end associate
  contains
    subroutine mapfull(buf, n, p)
      real(kind=fp), target, intent(inout) :: buf(:)
      integer, intent(in) :: n
      real(kind=fp), pointer, intent(out) :: p(:,:,:)
      p(1:n,1:n,1:3) => buf(1:n*n*3)
    end subroutine
  end subroutine

!-------------------------------------------------------------------------------
  subroutine giao_vxc(basis, molGrid, infos, coeffa, coeffb, urohf, &
                      vmata, vmatb, mxAngMom, nbf, dft_threshold)
    use mod_dft_molgrid, only: dft_grid_t
    use basis_tools, only: basis_set
    use mod_dft_gridint, only: xc_options_t, run_xc
    use types, only: information

    type(dft_grid_t), target, intent(in) :: molGrid
    type(information), target, intent(in) :: infos
    type(basis_set) :: basis
    logical, intent(in) :: urohf
    integer, intent(in) :: mxAngMom, nbf
    real(kind=fp), target, intent(inout) :: coeffa(nbf,*), coeffb(nbf,*)
    real(kind=fp), intent(out) :: vmata(3,nbf,nbf), vmatb(3,nbf,nbf)
    real(kind=fp), intent(in) :: dft_threshold

    type(xc_consumer_giao_t) :: dat
    type(xc_options_t) :: xc_opts
    integer :: i, j, t, ish, k, ao
    real(kind=fp), allocatable :: full(:,:,:)

    do j = 1, nbf
      coeffa(:,j) = coeffa(:,j)*basis%bfnrm(:)
    end do
    if (urohf) then
      do j = 1, nbf
        coeffb(:,j) = coeffb(:,j)*basis%bfnrm(:)
      end do
    end if

    allocate(dat%aocent(3,nbf))
    do ish = 1, basis%nshell
      do k = 1, basis%naos(ish)
        ao = basis%ao_offset(ish) + k - 1
        dat%aocent(:,ao) = basis%shell_centers(ish,1:3)
      end do
    end do

    xc_opts%isGGA = infos%functional%needGrd
    xc_opts%needTau = .false.
    xc_opts%functional => infos%functional
    xc_opts%hasBeta = urohf
    xc_opts%isWFVecs = .true.
    xc_opts%numAOs = nbf
    xc_opts%maxPts = molGrid%maxSlicePts
    xc_opts%limPts = molGrid%maxNRadTimesNAng
    xc_opts%numAtoms = infos%mol_prop%natom
    xc_opts%maxAngMom = mxAngMom
    xc_opts%nDer = 0
    xc_opts%numOccAlpha = infos%mol_prop%nelec_A
    xc_opts%numOccBeta = infos%mol_prop%nelec_B
    xc_opts%wfAlpha => coeffa(:,1:nbf)
    xc_opts%wfBeta => coeffb(:,1:nbf)
    xc_opts%molGrid => molGrid
    xc_opts%dft_threshold = dft_threshold
    xc_opts%ao_threshold = infos%dft%grid_ao_threshold
    xc_opts%ao_sparsity_ratio = infos%dft%grid_ao_sparsity_ratio
    if (infos%dft%grid_pruned) xc_opts%ao_sparsity_ratio = 0.0_fp

    call dat%pe%init(infos%mpiinfo%comm, infos%mpiinfo%usempi)
    call run_xc(xc_opts, dat, basis)

    do j = 1, nbf
      coeffa(:,j) = coeffa(:,j)/basis%bfnrm(:)
    end do
    if (urohf) then
      do j = 1, nbf
        coeffb(:,j) = coeffb(:,j)/basis%bfnrm(:)
      end do
    end if

    allocate(full(nbf,nbf,3))
    full = reshape(dat%va2(1:nbf*nbf*3,1), [nbf,nbf,3])
    do t = 1, 3
      do i = 1, nbf
        do j = 1, nbf
          vmata(t,i,j) = (full(i,j,t) - full(j,i,t))*basis%bfnrm(i)*basis%bfnrm(j)
        end do
      end do
    end do
    if (urohf) then
      full = reshape(dat%vb2(1:nbf*nbf*3,1), [nbf,nbf,3])
      do t = 1, 3
        do i = 1, nbf
          do j = 1, nbf
            vmatb(t,i,j) = (full(i,j,t) - full(j,i,t))*basis%bfnrm(i)*basis%bfnrm(j)
          end do
        end do
      end do
    end if
    deallocate(full)
    call dat%clean()
  end subroutine

end module mod_dft_gridint_giao
