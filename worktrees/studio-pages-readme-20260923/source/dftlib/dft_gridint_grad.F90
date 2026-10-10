module mod_dft_gridint_grad

  use precision, only: fp
  use mod_dft_gridint, only: xc_engine_t, xc_consumer_t
  use mod_dft_gridint, only: OQP_FUNTYP_LDA, OQP_FUNTYP_MGGA
  use mod_dft_gridint, only: compAtGradRho, compAtGradDRho, compAtGradTau

  implicit none

!-------------------------------------------------------------------------------

  type, extends(xc_consumer_t) :: xc_consumer_grad_t
    real(kind=fp), allocatable :: bfgrad(:,:,:)
    real(kind=fp), allocatable :: nucgrad(:,:,:)
    real(kind=fp), allocatable :: exc_weighted(:,:)
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
    real(kind=fp), allocatable :: tmp_(:,:)
    real(kind=fp), allocatable :: d1dsx(:,:,:) !< Temporary storage for dE/d\sigma
  contains
    procedure :: parallel_start
    procedure :: parallel_stop
    procedure :: resetGradPointers
    procedure :: update
    procedure :: postUpdate
    procedure :: clean
  end type

!-------------------------------------------------------------------------------

  private
  public derexc_blk

!-------------------------------------------------------------------------------

contains

!-------------------------------------------------------------------------------

  subroutine parallel_start(self, xce, nthreads)
    implicit none
    class(xc_consumer_grad_t), target, intent(inout) :: self
    class(xc_engine_t), intent(in) :: xce
    integer, intent(in) :: nthreads
    integer :: nat, i, j
    call clean_work(self)
    allocate( self%bfGrad(xce%numAOs, 3, nthreads) &
            , self%nucGrad(3, xce%numAtoms, nthreads) &
            , self%exc_weighted(xce%maxPts, nthreads) &
            , self%d1dsx(xce%maxPts, 3, nthreads) &
            , self%tmp_(xce%numAOs*3, nthreads) &
            , source=0.0d0)

    nat = xce%numAtoms
    allocate( self%part_rij(nat,nat) &
            , self%part_rhat(3,nat,nat) &
            , self%part_dist(nat,nthreads) &
            , self%part_cells(nat,nthreads) &
            , self%part_dlog(3,nat,nat,nthreads) &
            , self%part_dsum(3,nat,nthreads) &
            , source=0.0_fp)
    do i = 1, nat
      do j = 1, nat
        if (i == j) cycle
        self%part_rij(i,j) = norm2(self%atom_xyz(:,i)-self%atom_xyz(:,j))
        if (self%part_rij(i,j) > tiny(1.0_fp)) &
          self%part_rhat(:,i,j) = &
            (self%atom_xyz(:,i)-self%atom_xyz(:,j))/self%part_rij(i,j)
      end do
    end do
  end subroutine

!-------------------------------------------------------------------------------

  subroutine parallel_stop(self)
    implicit none
    class(xc_consumer_grad_t), intent(inout) :: self
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
    class(xc_consumer_grad_t), intent(inout) :: self
    if (allocated(self%bfGrad)) deallocate(self%bfGrad)
    if (allocated(self%nucGrad)) deallocate(self%nucGrad)
    if (allocated(self%exc_weighted)) deallocate(self%exc_weighted)
    if (allocated(self%d1dsx)) deallocate(self%d1dsx)
    if (allocated(self%tmp_)) deallocate(self%tmp_)
    if (allocated(self%part_rij)) deallocate(self%part_rij)
    if (allocated(self%part_rhat)) deallocate(self%part_rhat)
    if (allocated(self%part_dist)) deallocate(self%part_dist)
    if (allocated(self%part_cells)) deallocate(self%part_cells)
    if (allocated(self%part_dlog)) deallocate(self%part_dlog)
    if (allocated(self%part_dsum)) deallocate(self%part_dsum)
  end subroutine clean_work

  subroutine clean(self)
    implicit none
    class(xc_consumer_grad_t), intent(inout) :: self
    call clean_work(self)
    if (allocated(self%atom_xyz)) deallocate(self%atom_xyz)
    if (allocated(self%dummy_atom)) deallocate(self%dummy_atom)
    if (allocated(self%surface_shift)) deallocate(self%surface_shift)
  end subroutine clean

!-------------------------------------------------------------------------------
!> @brief Adjust internal memory storage for a given
!>  number of pruned grid points
!> @author Konstantin Komarov
 subroutine resetGradPointers(self, xce, tmp, myThread)
    class(xc_consumer_grad_t), target, intent(inout) :: self
    class(xc_engine_t), intent(in) :: xce
    real(kind=fp), intent(out), pointer :: tmp(:,:)
    integer, intent(in) :: myThread

!   pruned AOs or no pruned AOs
    associate ( numAOs => xce%numAOs_p &  ! number of pruned AOs
      )
      tmp(1:numAOs, 1:3) => self%tmp_(1:numAOs*3, myThread)
    end associate

 end subroutine

!-------------------------------------------------------------------------------

 subroutine update(self, xce, mythread)

    class(xc_consumer_grad_t), intent(inout) :: self
    class(xc_engine_t), intent(in) :: xce
    integer :: mythread, i
    real(kind=fp), pointer :: tmpGrad(:,:)

    call self%resetGradPointers(xce, tmpGrad,  myThread)

    associate ( bfgrad  => self%bfgrad(:,:,mythread) &
              , d1dsx   => self%d1dsx(:,:,mythread) &
              , aoG1    => xce%aoG1 &
              , aoG2    => xce%aoG2 &
              , moVA    => xce%moVA &
              , moVB    => xce%moVB &
              , moG1A   => xce%moG1A &
              , moG1B   => xce%moG1B &
              , hasBeta => xce%hasBeta &
              , numPts  => xce%numPts &
              , xc      => xce%XCLib &
              , drho    => xce%xclib%drho  &
              , ids     => xce%XCLib%ids &
              , d1ds    => xce%XCLib%d1ds &
              , d1dr    => xce%XCLib%d1dr &
              , d1dt    => xce%XCLib%d1dt &
      )
      tmpGrad = 0.0d0

!     LDA gradient
      call compAtGradRho(tmpGrad, d1dr(1,:), moVA, aoG1, numPts)

!     GGA gradient
      if (xce%funTyp /= OQP_FUNTYP_LDA) then
          do i = 1, numPts
            d1dsx(i,1:3) = 2*d1ds(ids%ga,i)*drho(1:3,i)+d1ds(ids%gc,i)*drho(4:6,i)
          end do

          call compAtGradDRho(tmpGrad, d1dsx, moVA, moG1A, aoG1, aoG2, numPts)
      end if

!     Meta-GGA gradient
      if (xce%funTyp == OQP_FUNTYP_MGGA) &
        call compAtGradTau(tmpGrad, d1dt(1,:), moG1A, aoG2, numPts)


      if (hasBeta) then
        call compAtGradRho(tmpGrad, d1dr(2,:), moVB, aoG1, numPts)

        if (xce%funTyp /= OQP_FUNTYP_LDA) then
            do i = 1, numPts
              d1dsx(i,1:3) = 2*d1ds(ids%gb,i)*drho(4:6,i)+d1ds(ids%gc,i)*drho(1:3,i)
            end do
            call compAtGradDRho(tmpGrad, d1dsx, moVB, moG1B, aoG1, aoG2, numPts)
        end if

        if (xce%funTyp == OQP_FUNTYP_MGGA) &
          call compAtGradTau(tmpGrad, d1dt(2,:), moG1B, aoG2, numPts)
      end if

      ! The native XC energy array is unweighted.  Retain the integrand of the
      ! finite molecular quadrature so its fuzzy-cell weight can be
      ! differentiated consistently with the energy evaluation.
      self%exc_weighted(1:numPts,mythread) = &
        xc%exc(1:numPts)*xce%xyzw(1:numPts,4)

      call add_partition_weight_gradient(self, xce, mythread)
      ! A radial/angular slice moves rigidly with its owner atom.  The AO
      ! spatial derivative in tmpGrad is the corresponding owner-point term.
      self%nucgrad(:,xce%currAtom,mythread) = &
        self%nucgrad(:,xce%currAtom,mythread) + sum(tmpGrad, dim=1)

    end associate
 end subroutine

!-------------------------------------------------------------------------------

!> Add the derivative of the normalized atom-centred fuzzy-cell weight to the
!> ground-state XC energy.  Together with the owner-point term accumulated in
!> update, this differentiates the same finite quadrature used for the energy.
 subroutine add_partition_weight_gradient(self, xce, mythread)
    use mod_dft_partfunc, only: partition_function

    class(xc_consumer_grad_t), intent(inout) :: self
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
    call partfunc%set(self%part_fun_type)

    associate ( dist => self%part_dist(:,mythread) &
              , cells => self%part_cells(:,mythread) &
              , dlog => self%part_dlog(:,:,:,mythread) &
              , dsum => self%part_dsum(:,:,mythread))
      do ipt = 1, xce%numPts
        q_weighted = self%exc_weighted(ipt,mythread)
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
            + q_weighted*(dlog(:,b,owner)-dfactor*dsum(:,b))
        end do
      end do
    end associate
 end subroutine add_partition_weight_gradient

 subroutine postUpdate(self, xce, mythread)

    class(xc_consumer_grad_t), intent(inout) :: self
    class(xc_engine_t), intent(in) :: xce
    integer :: mythread

    real(kind=fp), pointer :: tmpGrad(:,:)

    call self%resetGradPointers(xce, tmpGrad,  myThread)

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

!-------------------------------------------------------------------------------

!> @brief Compute grid XC contribution to the nuclear gradient
!> @note  Includes normalized partition-weight derivatives and rigid motion of
!>  each atom-centred radial/angular slice, so the analytic derivative matches
!>  the finite XC quadrature used for the energy.
!> @param[in]    da        density matrix, alpha-spin
!> @param[in]    db        density matrix, beta-spin
!> @param[inout] dedft     nuclear gradient
!> @param[out]   totele    electronic denisty integral
!> @param[out]   totkin    kinetic energy integral
!> @param[in]    mxAngMom  max. needed ang. mom. value (incl. derivatives)
!> @param[in]    nbf        basis set size
!> @param[in]    isGGA     .TRUE. if GGA/mGGA functional used
!> @param[in]    urohf     .TRUE. if open-shell calculation
!> @author Vladimir Mironov
  subroutine derexc_blk(basis, molGrid, da, db, dedft, &
                        totele, totkin, &
                        mxAngMom, nbf, dft_threshold, urohf, infos)
!$  use omp_lib, only: omp_get_num_threads, omp_get_thread_num
    use basis_tools, only: basis_set
    use mod_dft_gridint, only: xc_options_t, run_xc
    use types, only: information
    use mod_dft_molgrid, only: dft_grid_t

    implicit none

    type(information), target, intent(in) :: infos
    type(dft_grid_t), target, intent(in) :: molGrid

    type(basis_set) :: basis
    logical, intent(IN) :: urohf
    integer, intent(IN) :: mxAngMom, nbf
    real(KIND=fp), intent(INOUT) :: totele, totkin
    real(KIND=fp), intent(INOUT) :: da(nbf, *), db(nbf, *), dedft(:, :)
    real(kind=fp), intent(in) :: dft_threshold

    type(xc_consumer_grad_t) :: dat
    type(xc_options_t) :: xc_opts

    integer :: j

    integer :: nat

    real(KIND=fp), target, allocatable :: da2(:, :), db2(:, :)


    nat = infos%mol_prop%natom

    allocate (da2(nbf, nbf))
    do j = 1, nbf
      da2(:, j) = da(:, j)*basis%bfnrm(j)*basis%bfnrm(1:nbf)
    end do
    if (urohf) then
      allocate (db2(nbf, nbf))
      do j = 1, nbf
        db2(:, j) = db(:, j)*basis%bfnrm(j)*basis%bfnrm(1:nbf)
      end do
    end if

    xc_opts%isGGA = infos%functional%needGrd
    xc_opts%needTau = infos%functional%needTau
    xc_opts%functional => infos%functional
    xc_opts%hasBeta = urohf
    xc_opts%isWFVecs = .false.
    xc_opts%numAOs = nbf
    xc_opts%maxPts = molGrid%maxSlicePts
    xc_opts%limPts = molGrid%maxNRadTimesNAng
    xc_opts%numAtoms = infos%mol_prop%natom
    xc_opts%maxAngMom = mxAngMom
    xc_opts%nDer = 1
    xc_opts%numOccAlpha = infos%mol_prop%nelec_A
    xc_opts%numOccBeta = infos%mol_prop%nelec_B
    xc_opts%wfAlpha => da2
    xc_opts%wfBeta => db2
    xc_opts%dft_threshold = dft_threshold
    xc_opts%molGrid => molGrid

    dat%part_fun_type = molGrid%partFunType
    dat%has_surface_shift = molGrid%hasSurfaceShift
    dat%atom_xyz = infos%atoms%xyz
    dat%dummy_atom = molGrid%dummyAtom
    dat%surface_shift = molGrid%surfaceShift

    call dat%pe%init(infos%mpiinfo%comm, infos%mpiinfo%usempi)

    call run_xc(xc_opts, dat, basis)

    totele = dat%N_elec
    totkin = dat%E_kin

    do j = 1, basis%nshell
      associate (atom => basis%origin(j), &
                 offset => basis%ao_offset(j), &
                 naos => basis%naos(j))
        dedft(1, atom) = dedft(1, atom)-sum(dat%bfGrad(offset:offset+naos-1, 1, 1))
        dedft(2, atom) = dedft(2, atom)-sum(dat%bfGrad(offset:offset+naos-1, 2, 1))
        dedft(3, atom) = dedft(3, atom)-sum(dat%bfGrad(offset:offset+naos-1, 3, 1))
      end associate
    end do

    dedft = dedft + dat%nucGrad(:,:,1)

    deallocate (da2)
    if (urohf) deallocate (db2)

    call dat%clean()

  end subroutine

!-------------------------------------------------------------------------------

end module mod_dft_gridint_grad
