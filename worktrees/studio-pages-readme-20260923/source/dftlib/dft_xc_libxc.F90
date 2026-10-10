 module mod_dft_xc_libxc
    use precision, only: fp
    use mod_dft_xclib
    implicit none

    type, extends(xc_lib_t) :: xc_libxc_t
        integer :: nSpin
        real(kind=fp), contiguous, pointer :: lib_output(:)  => null()
    contains
        procedure init
        procedure setPts
        procedure compute
    end type

 contains

 subroutine init(self, reqSigma, reqTau, reqLapl, reqBeta, maxPts, nDer)
    ! E_XC and E_CORR are separate
    class(xc_libxc_t) :: self
    logical, intent(in) :: reqSigma, reqTau, reqLapl, reqBeta
    integer, intent(in) :: maxPts, nDer

    integer :: ndens

    call self%clean()

    self%reqSigma = reqSigma
    self%reqTau   = reqTau
    self%reqLapl  = reqLapl
    self%reqBeta  = reqBeta
    self%maxPts   = maxPts
    self%nDer     = nDer

    self%providesEXC  = .TRUE.
    self%providesEX   = .FALSE.
    self%providesEC   = .FALSE.

    self%nSpin = 1
    if (self%reqBeta) self%nSpin = 2

    ndens = 2 + 6 + 3 + 2 + 2 &
      + 1 + 2 + 3 + 2 + 2

    if (nDer > 1) &
      ndens = ndens &
        + 3 + 6 + 4 + 4 + 6 + 6 + 6 + 3 + 4 + 3

    if (nDer > 2) &
      ndens = ndens &
        + 4 + 10 + 9 + 12 + 4 + 6 + 6 + 12 + 12 &
        + 9 + 6 + 6 + 12 + 8 + 12 + 9 + 12 + 6 &
        + 6 + 4

    allocate(self%memory_(1:maxPts*ndens))

    call self%resetEnergy

 end subroutine

 subroutine setPts(self, numPts)
    class(xc_libxc_t), target :: self
    integer, intent(in) :: numPts
    integer :: i, i_out0


    self%numPts = numPts

    i = 0

    ! Input
    self%rho  => addmem(self%memory_, i, 2,   numPts)
    self%drho => addmem(self%memory_, i, 6,   numPts)

    self%sig  => addmem(self%memory_, i, 3,   numPts)
    self%tau  => addmem(self%memory_, i, 2,   numPts)
    self%lapl => addmem(self%memory_, i, 2,   numPts)

    ! Output
    i_out0 = i+1

    self%exc(1:numPts) => self%memory_(i+1:)
    i = i + numPts

    self%d1dr => addmem(self%memory_, i, 2,   numPts)
    self%d1ds => addmem(self%memory_, i, 3,   numPts)
    self%d1dt => addmem(self%memory_, i, 2,   numPts)
    self%d1dl => addmem(self%memory_, i, 2,   numPts)

    self%lib_output => self%memory_(i_out0:i)
    if (self%nDer==1) return

    self%d2r2 => addmem(self%memory_, i, 3,   numPts)
    self%d2rs => addmem(self%memory_, i, 6,   numPts)
    self%d2rt => addmem(self%memory_, i, 4,   numPts)
    self%d2rl => addmem(self%memory_, i, 4,   numPts)
    self%d2s2 => addmem(self%memory_, i, 6,   numPts)
    self%d2st => addmem(self%memory_, i, 6,   numPts)
    self%d2sl => addmem(self%memory_, i, 6,   numPts)
    self%d2t2 => addmem(self%memory_, i, 3,   numPts)
    self%d2tl => addmem(self%memory_, i, 4,   numPts)
    self%d2l2 => addmem(self%memory_, i, 3,   numPts)

    self%lib_output => self%memory_(i_out0:i)
    if (self%nDer==2) return

    self%d3r3  => addmem(self%memory_, i, 4,  numPts)
    self%d3s3  => addmem(self%memory_, i, 10, numPts)
    self%d3r2s => addmem(self%memory_, i, 9,  numPts)
    self%d3rs2 => addmem(self%memory_, i, 12, numPts)

    self%d3t3  => addmem(self%memory_, i, 4,  numPts)
    self%d3r2t => addmem(self%memory_, i, 6,  numPts)
    self%d3rt2 => addmem(self%memory_, i, 6,  numPts)
    self%d3rst => addmem(self%memory_, i, 12, numPts)
    self%d3s2t => addmem(self%memory_, i, 12, numPts)
    self%d3st2 => addmem(self%memory_, i, 9,  numPts)

    self%d3r2l => addmem(self%memory_, i, 6,  numPts)
    self%d3rl2 => addmem(self%memory_, i, 6,  numPts)
    self%d3rsl => addmem(self%memory_, i, 12, numPts)
    self%d3rtl => addmem(self%memory_, i, 8,  numPts)
    self%d3s2l => addmem(self%memory_, i, 12, numPts)
    self%d3sl2 => addmem(self%memory_, i, 9,  numPts)
    self%d3stl => addmem(self%memory_, i, 12, numPts)
    self%d3t2l => addmem(self%memory_, i, 6,  numPts)
    self%d3tl2 => addmem(self%memory_, i, 6,  numPts)
    self%d3l3  => addmem(self%memory_, i, 4,  numPts)

    self%lib_output  => self%memory_(i_out0:i)

 contains

   function addmem(memory, pos, d1, d2) result(res)
     real(kind=fp), contiguous, target :: memory(:)
     integer, intent(inout) :: pos
     integer, intent(in) :: d1, d2
     real(kind=fp), contiguous, pointer :: res(:,:)
     res(1:d1,1:d2) => memory(pos+1: )
     pos = pos + d1*d2
   end function

 end subroutine

 subroutine compute(self, functional, wts)
    use functionals, only: functional_t
!    use xc_f03_lib_m
    class(xc_libxc_t) :: self
    class(functional_t) :: functional
    real(kind=fp), intent(in) :: wts(:)

    self%lib_output = 0

    select case (self%nDer)
    case (1)
      call functional%calc_evxc(self%numPts, &
              rho      = self%rho, &
              sigma    = self%sig, &
              tau      = self%tau, &
              lapl     = self%lapl, &
              energy   = self%exc, &
              dedrho   = self%d1dr, &
              dedsigma = self%d1ds, &
              dedtau   = self%d1dt, &
              dedlapl  = self%d1dl)
    case (2)
      call functional%calc_evfxc(self%numPts, &
            rho         = self%rho, &
            sigma       = self%sig, &
            tau         = self%tau, &
            lapl        = self%lapl, &
            energy      = self%exc, &
            dedrho      = self%d1dr, &
            dedsigma    = self%d1ds, &
            dedtau      = self%d1dt, &
            dedlapl     = self%d1dl, &
            v2rho2      = self%d2r2, &
            v2sigma2    = self%d2s2, &
            v2tau2      = self%d2t2, &
            v2lapl2     = self%d2l2, &
            v2rhosigma  = self%d2rs, &
            v2rhotau    = self%d2rt, &
            v2rholapl   = self%d2rl, &
            v2sigmatau  = self%d2st, &
            v2sigmalapl = self%d2sl, &
            v2lapltau   = self%d2tl)
    case (3)
      call functional%calc_xc(self%numPts, &
            rho            = self%rho, &
            sigma          = self%sig, &
            tau            = self%tau, &
            lapl           = self%lapl, &
            energy         = self%exc, &
            dedrho         = self%d1dr, &
            dedsigma       = self%d1ds, &
            dedlapl        = self%d1dl, &
            dedtau         = self%d1dt, &
            v2rho2         = self%d2r2, &
            v2rhosigma     = self%d2rs, &
            v2rholapl      = self%d2rl, &
            v2rhotau       = self%d2rt, &
            v2sigma2       = self%d2s2, &
            v2sigmalapl    = self%d2sl, &
            v2sigmatau     = self%d2st, &
            v2lapl2        = self%d2l2, &
            v2lapltau      = self%d2tl, &
            v2tau2         = self%d2t2, &
            v3rho3         = self%d3r3, &
            v3rho2sigma    = self%d3r2s, &
            v3rho2lapl     = self%d3r2l, &
            v3rho2tau      = self%d3r2t, &
            v3rhosigma2    = self%d3rs2, &
            v3rhosigmalapl = self%d3rsl, &
            v3rhosigmatau  = self%d3rst, &
            v3rholapl2     = self%d3rl2, &
            v3rholapltau   = self%d3rtl, &
            v3rhotau2      = self%d3rt2, &
            v3sigma3       = self%d3s3, &
            v3sigma2lapl   = self%d3s2l, &
            v3sigma2tau    = self%d3s2t, &
            v3sigmalapl2   = self%d3sl2, &
            v3sigmalapltau = self%d3stl, &
            v3sigmatau2    = self%d3st2, &
            v3lapl3        = self%d3l3, &
            v3lapl2tau     = self%d3tl2, &
            v3lapltau2     = self%d3t2l, &
            v3tau3         = self%d3t3)
    end select

    self%E_xc   = self%E_xc   + dot_product(self%exc, wts)

    call self%scaleXC(wts)


 end subroutine

 end module
