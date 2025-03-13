module int2_pairs

  use precision, only: dp
  implicit none

  private
  public int2_cutoffs_t
  public int2_pair_storage

  real(dp), parameter :: pi4 = atan(1.0_dp) !0.78539816339744831_dp
  real(dp), parameter :: pi = 4*pi4
  real(dp), parameter :: sqrtpito52 = sqrt(2.0d0)*pi**1.25d0

  type int2_pair_storage
      real(kind=dp), allocatable :: &
        alpha_a(:), alpha_b(:), g(:), ginv(:), &
        k(:), p(:,:), pa(:,:), pb(:,:)
      real(kind=dp), allocatable :: &
        rab(:), uab(:)
      integer, allocatable :: ppid(:,:)
  contains
      procedure :: alloc => int2_prepare_pair_storage
      procedure :: compute => int2_prepare_shellpairs
      procedure :: clean   => int2_clean_pair_storage
  end type

  type int2_cutoffs_t
     real(dp) :: integral_cutoff
     real(dp) :: pair_cutoff, quartet_cutoff, exponent_cutoff
     real(dp) :: pair_cutoff_squared
     real(dp) :: quartet_cutoff_squared
  contains
    procedure :: get => get_int2_accuracy
    procedure :: set => set_int2_accuracy
  end type int2_cutoffs_t

contains

  subroutine int2_prepare_shellpairs(ppairs, basis, cutoffs, noswap)
    use basis_tools, only: basis_set
    implicit none
    class(int2_pair_storage), intent(inout) :: ppairs
    type(basis_set), intent(in) :: basis
    type(int2_cutoffs_t), intent(in) :: cutoffs
    logical, optional, intent(in) :: noswap
    integer :: i, j

    do i = 1, basis%nshell
        do j = 1, i
            call int2_prepare_pair(ppairs, basis, cutoffs, i, j, noswap)
        end do
    end do
  end subroutine

! Cleanup
  subroutine int2_clean_pair_storage(ppairs)
    implicit none
    class(int2_pair_storage), intent(inout) :: ppairs
    if (allocated(ppairs%alpha_a)) deallocate(ppairs%alpha_a)
    if (allocated(ppairs%alpha_b)) deallocate(ppairs%alpha_b)
    if (allocated(ppairs%g)      ) deallocate(ppairs%g)
    if (allocated(ppairs%ginv)   ) deallocate(ppairs%ginv)
    if (allocated(ppairs%k)      ) deallocate(ppairs%k)
    if (allocated(ppairs%p)      ) deallocate(ppairs%p)
    if (allocated(ppairs%pa)     ) deallocate(ppairs%pa)
    if (allocated(ppairs%pb)     ) deallocate(ppairs%pb)
    if (allocated(ppairs%rab)    ) deallocate(ppairs%rab)
    if (allocated(ppairs%uab)    ) deallocate(ppairs%uab)
    if (allocated(ppairs%ppid)   ) deallocate(ppairs%ppid)
  end subroutine


! Count the required storage size
  subroutine int2_prepare_pair_storage(ppairs, basis, cutoffs)
    use basis_tools, only: basis_set
    implicit none
    class(int2_pair_storage), intent(inout) :: ppairs
    type(basis_set), intent(in) :: basis
    type(int2_cutoffs_t), intent(in) :: cutoffs

    integer :: sha, shb, ksa, ksb, ncona, nconb
    real(kind=dp) :: a(3), b(3), ab2
    integer :: p1, p2, p12
    real(kind=dp) :: alpha1_, alpha2_, gammap_, gpinv_, e12, klog
    real(kind=dp) :: logcut
    real(kind=dp), allocatable, target :: cclog(:)
    real(kind=dp), pointer :: cclogpa(:), cclogpb(:)

    integer :: nprim, nshell, nshell2
    integer :: pair_id
    integer :: npptotal

    nshell = basis%nshell
    nshell2 = nshell*(nshell+1)/2
    nprim = basis%g_offset(nshell) + basis%ncontr(nshell) - 1

    if (.not.allocated(ppairs%ppid).or.ubound(ppairs%ppid,2) < nshell2) then
      if (allocated(ppairs%ppid)) deallocate(ppairs%ppid)
      allocate(ppairs%ppid(2,nshell2))
    end if

!   Allocate storage for screening
    allocate(cclog(nprim))
    cclog = log(abs(basis%cc(:nprim)) )
    logcut = log(cutoffs%quartet_cutoff)

    npptotal = 0

    pair_id = 0
    do sha = 1, basis%nshell
      do shb = 1, sha

        ksa = basis%g_offset(sha)
        ksb = basis%g_offset(shb)
        ncona = basis%ncontr(sha)
        nconb = basis%ncontr(shb)

        cclogpa => cclog(ksa:ksa+ncona-1)
        cclogpb => cclog(ksb:ksb+nconb-1)

        a = basis%shell_centers(sha,1:3)
        b = basis%shell_centers(shb,1:3)
        AB2 = SUM((A - B)*(A - B))

        p12 = 0
        do p1 = 1, ncona
          do p2 = 1, nconb
            alpha1_ = basis%ex(ksa-1+p1)
            alpha2_ = basis%ex(ksb-1+p2)
            gammap_ = alpha1_ + alpha2_
            e12 = alpha1_*alpha2_*AB2
            if (e12 > gammap_*cutoffs%exponent_cutoff) cycle
            gpinv_ = 1/gammap_
            e12 = e12*gpinv_
            klog = cclogpa(p1) + cclogpb(p2) - e12
            if (klog < logcut) cycle
            p12 = p12 + 1
          end do
        end do

        pair_id = pair_id + 1
        ppairs%ppid(1, pair_id) = p12
        ppairs%ppid(2, pair_id) = npptotal + 1
        npptotal = npptotal + p12

      end do
    end do

    deallocate(cclog)

!   Allocate storage
    if (allocated(ppairs%alpha_a)) deallocate(ppairs%alpha_a)
    if (allocated(ppairs%alpha_b)) deallocate(ppairs%alpha_b)
    if (allocated(ppairs%g      )) deallocate(ppairs%g      )
    if (allocated(ppairs%ginv   )) deallocate(ppairs%ginv   )
    if (allocated(ppairs%k      )) deallocate(ppairs%k      )
    if (allocated(ppairs%p      )) deallocate(ppairs%p      )
    if (allocated(ppairs%pa     )) deallocate(ppairs%pa     )
    if (allocated(ppairs%pb     )) deallocate(ppairs%pb     )
    if (allocated(ppairs%rab    )) deallocate(ppairs%rab    )
    if (allocated(ppairs%uab    )) deallocate(ppairs%uab    )

    allocate(ppairs%alpha_a(  npptotal), source=0.0d0)
    allocate(ppairs%alpha_b(  npptotal), source=0.0d0)
    allocate(ppairs%g      (  npptotal), source=0.0d0)
    allocate(ppairs%ginv   (  npptotal), source=0.0d0)
    allocate(ppairs%k      (  npptotal), source=0.0d0)
    allocate(ppairs%p      (3,npptotal), source=0.0d0)
    allocate(ppairs%pa     (3,npptotal), source=0.0d0)
    allocate(ppairs%pb     (3,npptotal), source=0.0d0)
    allocate(ppairs%rab    (  npptotal), source=0.0d0)
    allocate(ppairs%uab    (  npptotal), source=0.0d0)

!    write(*,*) 'npptotal=', npptotal

  end subroutine


  subroutine int2_prepare_pair(ppairs, basis, cutoffs, sh1, sh2, noswap)
    use basis_tools, only: basis_set
    implicit none
    class(int2_pair_storage), intent(inout) :: ppairs
    type(basis_set), intent(in) :: basis
    type(int2_cutoffs_t), intent(in) :: cutoffs
    integer, intent(in) :: sh1, sh2
    logical, optional, intent(in) :: noswap

    integer :: ksa, ksb, ncona, nconb, anga, angb, sha, shb
    integer :: anga_, angb_
    real(kind=dp) :: a(3), b(3), ab2
    integer :: p1, p2, p12
    integer :: id1, id2, id12
    real(kind=dp) :: alpha1_, alpha2_, gammap_, gpinv_, e12, k1_
    logical :: noswap_

    anga_ = basis%am(sh1)
    angb_ = basis%am(sh2)

    noswap_ = .false.
    if (present(noswap)) noswap_ = noswap

    if (.not.noswap_ .and. anga_>angb_) then
      sha = sh2
      shb = sh1
      anga = angb_
      angb = anga_
    else
      sha = sh1
      shb = sh2
      anga = anga_
      angb = angb_
    end if


    ksa = basis%g_offset(sha)
    ksb = basis%g_offset(shb)
    ncona = basis%ncontr(sha)
    nconb = basis%ncontr(shb)

    a = basis%shell_centers(sha,1:3)
    b = basis%shell_centers(shb,1:3)


    AB2 = SUM((A - B)*(A - B))

    id1 = max(sha, shb)
    id2 = min(sha, shb)
    id12 = id1*(id1-1)/2 + id2
    p12 = ppairs%ppid(2, id12)

    if (ppairs%ppid(1, id12) > 0) then
      ppairs%rab(p12) = sqrt(AB2)
      ppairs%uab(p12) = 1.0d0/sqrt(AB2)
    end if

    associate( c1 => basis%cc(ksa:ksa+ncona-1) &
             , c2 => basis%cc(ksb:ksb+nconb-1) &
       )
      DO p1 = 1, ncona
         DO p2 = 1, nconb
           alpha1_ = basis%ex(ksa-1+p1)
           alpha2_ = basis%ex(ksb-1+p2)
           gammap_ = alpha1_ + alpha2_
           e12 = alpha1_*alpha2_*AB2
           if (e12 > gammap_*cutoffs%exponent_cutoff) cycle
           gpinv_ = 1/gammap_
           e12 = e12*gpinv_
           K1_ = c1(p1)*c2(p2)*exp(-e12)
           if (abs(k1_) < cutoffs%quartet_cutoff) cycle

           ppairs%alpha_a(p12) = alpha1_
           ppairs%alpha_b(p12) = alpha2_
           ppairs%g(p12) = gammap_
           ppairs%ginv(p12) = gpinv_

           ppairs%P(:,p12) = (alpha1_*A + alpha2_*B)*gpinv_
           ppairs%PA(:,p12) = ppairs%P(:,p12) - A
           ppairs%PB(:,p12) = ppairs%P(:,p12) - B
           ppairs%K(p12) = sqrtpito52*k1_

           p12 = p12 + 1
         END DO
      END DO
    end associate

  end subroutine

!> @brief Get screening parameters for rotated axis integral code:
!> - prefactor for single bra/ket shell pair;
!> - total prefactor;
!> - value of exponential coefficient.
  subroutine get_int2_accuracy(this, cutoff_integral_value, cutoff_prefactor_p, cutoff_prefactor_pq, cutoff_exp)
    class(int2_cutoffs_t), intent(in) :: this
    real(kind=dp), intent(out) :: cutoff_integral_value, cutoff_prefactor_p, cutoff_prefactor_pq, cutoff_exp
    cutoff_integral_value = this%integral_cutoff
    cutoff_prefactor_p  = this%pair_cutoff
    cutoff_prefactor_pq = this%quartet_cutoff
    cutoff_exp          = this%exponent_cutoff
  end subroutine

!> @brief Set screening parameters for rotated axis integral code:
!> - prefactor for single bra/ket shell pair;
!> - total prefactor;
!> - value of exponential coefficient.
  subroutine set_int2_accuracy(this, cutoff_integral_value, cutoff_prefactor_p, cutoff_prefactor_pq, cutoff_exp)
    class(int2_cutoffs_t), intent(inout) :: this
    real(kind=dp), intent(in) :: cutoff_integral_value, cutoff_prefactor_p, cutoff_prefactor_pq, cutoff_exp
    this%integral_cutoff        = cutoff_integral_value
    this%pair_cutoff            = cutoff_prefactor_p
    this%pair_cutoff_squared    = cutoff_prefactor_p**2
    this%quartet_cutoff         = cutoff_prefactor_pq
    this%quartet_cutoff_squared = cutoff_prefactor_pq**2
    this%exponent_cutoff        = cutoff_exp
  end subroutine

end module
