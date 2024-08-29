 module mod_dft_xclib

    use precision, only: fp
    use functionals, only: functional_t
    implicit none

    private

    integer, public, parameter :: XCLIB_LIBXC   = 0 ! External LIBXC library

    integer, parameter  :: NDENS_XC  = 18+10
    ! see tddfun, up to 2nd derivatives nxdim(2)=18, ncdim(2)=35
    ! E_XC and E_CORR are separate
    ! two additional arrays are needed for EX0 and EC0
    integer, parameter  :: NDENS_TD  = 18+18+35+2
    ! LibXC, same, but E_XC and E_CORR are summed up
    integer, parameter  :: NDENS_LXC = 18+35

    real(kind=fp), parameter :: &
        ZERO = 0.0D+00, TWO = 2.0D+00, HALF = 0.5D+00

    ! indices of xc arrays
    type, public :: xc_pack_t

        integer :: &
            ra = 1, rb = 2, &

            ga = 1, gc = 2, gb = 3, &

            ta = 1, tb = 2, &

            rara = 1, rarb = 2, rbrb = 3, &

            raga = 1, ragc = 2, ragb = 3, rbga = 4, rbgc = 5, rbgb = 6, &

            rata = 1, ratb = 2, rbta = 3, rbtb = 4, &

            gaga = 1, gagc = 2, gagb = 3, gcgc = 4, gbgc = 5, gbgb = 6, &

            gata = 1, gatb = 2, gcta = 3, gctb = 4, gbta = 5, gbtb = 6, &

            tata = 1, tatb = 2, tbtb = 3, &

            rarara = 1, rararb = 2, rarbrb = 3, rbrbrb = 4, &

            gagaga = 1, gagagc = 2, gagagb = 3, gagcgc = 4, gagbgc = 5, &
            gagbgb = 6, gcgcgc = 7, gbgcgc = 8, gbgbgc = 9, gbgbgb = 10, &

            raraga = 1, raragc = 2, raragb = 3, rarbga = 4, rarbgc = 5, &
            rarbgb = 6, rbrbga = 7, rbrbgc = 8, rbrbgb = 9, &

            ragaga = 1, ragagc = 2, ragagb = 3, ragcgc = 4, ragbgc = 5, ragbgb = 6, &
            rbgaga = 7, rbgagc = 8, rbgagb = 9, rbgcgc = 10, rbgbgc = 11, rbgbgb = 12, &

            tatata = 1, tatatb = 2, tatbtb = 3, tbtbtb = 4, &

            rarata = 1, raratb = 2, rarbta = 3, rarbtb = 4, rbrbta = 5, rbrbtb = 6, &

            ratata = 1, ratatb = 2, ratbtb = 3, rbtata = 4, rbtatb = 5, rbtbtb = 6, &

            ragata = 1, ragatb = 2, ragcta = 3, ragctb = 4, ragbta = 5, ragbtb = 6, &
            rbgata = 7, rbgatb = 8, rbgcta = 9, rbgctb = 10, rbgbta = 11, rbgbtb = 12, &

            gagata = 1, gagatb = 2, gagcta = 3, gagctb = 4, gagbta = 5, gagbtb = 6, &
            gcgcta = 7, gcgctb = 8, gbgcta = 9, gbgctb = 10, gbgbta = 11, gbgbtb = 12, &

            gatata = 1, gatatb = 2, gatbtb = 3, gctata = 4, gctatb = 5, gctbtb = 6, &
            gbtata = 7, gbtatb = 8, gbtbtb = 9

    end type

    type, abstract, public :: xc_lib_t

        logical :: reqSigma  = .FALSE.  !< needs density gradients
        logical :: reqTau    = .FALSE.  !< needs k.e. density gradients
        logical :: reqLapl   = .FALSE.  !< needs laplacian of the density
        logical :: reqBeta   = .FALSE.  !< UHF flag

        integer :: maxPts    = 0
        integer :: numPts    = 0
        integer :: nDer      = 0
!        integer :: funTyp    = 0   !< 0 - LDA, 1 - GGA, 2 - MGGA

        real(kind=fp) :: E_xc       = 0.0
        real(kind=fp) :: E_exch     = 0.0
        real(kind=fp) :: E_corr     = 0.0

        logical :: providesEXC   = .FALSE.  !< Can get E_xc
        logical :: providesEX    = .FALSE.  !< Can get E_exch
        logical :: providesEC    = .FALSE.  !< Can get E_corr

!       Library id
        integer :: xclibID = XCLIB_LIBXC

        type(xc_pack_t) :: ids

        real(kind=fp), allocatable :: memory_(:)

        real(kind=fp), contiguous, pointer :: &
!           Input data
              rho(:,:)      => NULL() & !< density
            , drho(:,:)     => NULL() & !< gradient density
            , sig(:,:)      => NULL() & !< contracted density gradient
            , tau(:,:)      => NULL() & !< K.E. density
            , lapl(:,:)     => NULL() & !< Laplacian of the density

!           Output data
            , exc(:)        => NULL() & !< E(XC)
            , d1dr(:,:)     => NULL() & !< E(XC) LDA values
            , d1ds(:,:)     => NULL() & !< E(XC) GGA values
            , d1dt(:,:)     => NULL() & !< E(XC) MGGA values
            , d1dl(:,:)     => NULL() & !< E(XC) Laplacian

            , d2r2(:,:)     => NULL() & !< second derivatives of functional vs rho^2
            , d2s2(:,:)     => NULL() & !< second derivatives of functional vs sigma^2
            , d2t2(:,:)     => NULL() & !< second derivatives of functional vs tau^2
            , d2rs(:,:)     => NULL() & !< second derivatives of functional vs rho and sigma
            , d2rt(:,:)     => NULL() & !< second derivatives of functional vs rho and tau
            , d2st(:,:)     => NULL() & !< second derivatives of functional vs sigma and tau

            , d2rl(:,:)     => NULL() & !< second derivatives of functional vs rho and lapl
            , d2sl(:,:)     => NULL() & !< second derivatives of functional vs sigma and lapl
            , d2tl(:,:)     => NULL() & !< second derivatives of functional vs tau and lapl
            , d2l2(:,:)     => NULL() & !< second derivatives of functional vs lapl^2

            , d3r3(:,:)     => NULL() & !< third derivatives of functional vs rho^3
            , d3r2s(:,:)    => NULL() & !< third derivatives of functional vs rho^2 and sigma
            , d3rs2(:,:)    => NULL() & !< third derivatives of functional vs rho and sigma^2
            , d3s3(:,:)     => NULL() & !< third derivatives of functional vs sigma^3
            , d3t3(:,:)     => NULL() & !< third derivatives of functional vs tau^3
            , d3r2t(:,:)    => NULL() & !< third derivatives of functional vs rho^2 and tau
            , d3s2t(:,:)    => NULL() & !< third derivatives of functional vs sigma^2 and tau
            , d3rt2(:,:)    => NULL() & !< third derivatives of functional vs rho and tau^2
            , d3st2(:,:)    => NULL() & !< third derivatives of functional vs sigma and tau^2
            , d3rst(:,:)    => NULL() & !< third derivatives of functional vs rho, sigma and tau
            , d3r2l(:,:)    => NULL() & !<
            , d3rl2(:,:)    => NULL() & !<
            , d3rsl(:,:)    => NULL() & !<
            , d3rtl(:,:)    => NULL() & !<
            , d3s2l(:,:)    => NULL() & !<
            , d3sl2(:,:)    => NULL() & !<
            , d3stl(:,:)    => NULL() & !<
            , d3t2l(:,:)    => NULL() & !<
            , d3tl2(:,:)    => NULL() & !<
            , d3l3(:,:)     => NULL()

    contains
        procedure(init_xc_lib),     deferred :: init
        procedure(compute_xc_lib),  deferred :: compute
        procedure(setPts_xc_lib),   deferred :: setPts
        procedure :: clean
        procedure :: scalexc
        procedure, non_overridable :: echo
        procedure, non_overridable :: getEnergy
        procedure, non_overridable :: resetEnergy
    end type xc_lib_t

    abstract interface
        subroutine init_xc_lib(self, reqSigma, reqTau, reqLapl, reqBeta, maxPts, nDer)
            import
            class(xc_lib_t) :: self
            logical, intent(in) :: reqSigma, reqTau, reqLapl, reqBeta
            integer, intent(in) :: maxPts, nDer
        end subroutine
        subroutine setPts_xc_lib(self, numPts)
            import
            class(xc_lib_t), target :: self
            integer, intent(in) :: numPts
        end subroutine
        subroutine compute_xc_lib(self, functional, wts)
            import
            class(xc_lib_t) :: self
            class(functional_t) :: functional
            real(kind=fp), intent(in) :: wts(:)
        end subroutine
    end interface

 contains

!> @brief Print parameters of the xc_engine_t instance
!> @author Vladimir Mironov
 subroutine echo(self)
    class(xc_lib_t) :: self
    write(*,*) 'reqSigma =', self%reqSigma
    write(*,*) 'reqTau   =', self%reqTau
    write(*,*) 'reqBeta  =', self%reqBeta
    write(*,*) 'maxPts   =', self%maxPts
    write(*,*) 'numPts   =', self%numPts
    write(*,*) 'nDer     =', self%nDer
    write(*,*) 'xclibID  =', self%xclibID

 end subroutine

!> @brief Get debug statistics
!> @author Vladimir Mironov
 subroutine getEnergy(self, E_xc, E_exch, E_corr)
    class(xc_lib_t) :: self
    real(kind=fp), intent(out) :: &
        E_xc, E_exch, E_corr

    E_xc       = self%E_xc
    E_exch     = self%E_exch
    E_corr     = self%E_corr

 end subroutine

!> @brief Set debug statistics
!> @author Vladimir Mironov
 subroutine resetEnergy(self)
    class(xc_lib_t) :: self

    self%E_exch     = 0.0
    self%E_corr     = 0.0

 end subroutine

!> @brief Cleanup
!> @author Vladimir Mironov
 subroutine clean(self)
    class(xc_lib_t) :: self

    if (allocated(self%memory_)) deallocate(self%memory_)

    self%rho      => NULL()
    self%sig      => NULL()
    self%tau      => NULL()
    self%lapl     => NULL()
    self%exc      => NULL()
    self%d1dr     => NULL()
    self%d1ds     => NULL()
    self%d1dt     => NULL()
    self%d1dl     => NULL()
    self%d2r2     => NULL()
    self%d2s2     => NULL()
    self%d2t2     => NULL()
    self%d2rs     => NULL()
    self%d2rt     => NULL()
    self%d2st     => NULL()
    self%d2rl     => NULL()
    self%d2sl     => NULL()
    self%d2tl     => NULL()
    self%d2l2     => NULL()
    self%d3r3     => NULL()
    self%d3r2s    => NULL()
    self%d3rs2    => NULL()
    self%d3s3     => NULL()
    self%d3t3     => NULL()
    self%d3r2t    => NULL()
    self%d3s2t    => NULL()
    self%d3rt2    => NULL()
    self%d3st2    => NULL()
    self%d3rst    => NULL()
    self%d3r2l    => NULL()
    self%d3rl2    => NULL()
    self%d3rsl    => NULL()
    self%d3rtl    => NULL()
    self%d3s2l    => NULL()
    self%d3sl2    => NULL()
    self%d3stl    => NULL()
    self%d3t2l    => NULL()
    self%d3tl2    => NULL()
    self%d3l3     => NULL()

 end subroutine

!> @brief Scale XC values by grid weights
!> @author Vladimir Mironov
 subroutine scalexc(self, wts)
    class(xc_lib_t) :: self
    real(kind=fp) :: wts(:)

    call scale_2d(self%d1dr, wts)

    if (self%reqSigma) then
        call scale_2d(self%d1ds, wts)
    end if

    if (self%reqTau) then
        call scale_2d(self%d1dt, wts)
    end if

    if (self%nDer<2) return

    call scale_2d(self%d2r2, wts)

    if (self%reqSigma) then
        call scale_2d(self%d2s2, wts)
        call scale_2d(self%d2rs, wts)
    end if

    if (self%reqTau) then
        call scale_2d(self%d2t2, wts)
        call scale_2d(self%d2rt, wts)
        call scale_2d(self%d2st, wts)
    end if

    if (self%nDer<3) return

    call scale_2d(self%d3r3, wts)

    if (self%reqSigma) then
      call scale_2d(self%d3r2s, wts)
      call scale_2d(self%d3rs2, wts)
      call scale_2d(self%d3s3, wts)
    end if

    if (self%reqTau) then
      call scale_2d(self%d3r2t, wts)
      call scale_2d(self%d3rst, wts)
      call scale_2d(self%d3rt2, wts)
      call scale_2d(self%d3s2t, wts)
      call scale_2d(self%d3st2, wts)
      call scale_2d(self%d3t3, wts)
    end if

 end subroutine

!> @brief Scale 2d array along 1st dimension by a given
!>  vector of weights
 subroutine scale_2d(array, weights)
    real(kind=fp), intent(inout) :: array(:,:)
    real(kind=fp), intent(in) :: weights(:)
    integer :: i
    do i = lbound(array,1), ubound(array,1)
        array(i,:) = array(i,:) * weights
    end do
 end subroutine

 end module
