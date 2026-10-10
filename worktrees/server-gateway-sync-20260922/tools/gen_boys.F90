module boys_prepare
  implicit none
  integer, parameter :: dp = 8
  integer, parameter :: qp = 16
  integer, parameter :: &
    ntx = 4, &    !< max order of interpolation
    npf = 450, &  !< number of grid points for gamma
    ngrd = 7, &   !< number of grids for gamma
    npx = 1000, & !< number of grid points for exp(-x)
    mxqt = 16     !< max l of quartet

  integer, parameter :: &
    ntx_al = 7    !< max order of interpolation

  integer, parameter :: &
    igrid(0:16) = (/0, 1, 2, 3, 4, 5, 5, 5, 5, 6, 6, 6, 6, 7, 7, 7, 7/), &
    irgrd(0:16) = (/0, 1, 2, 3, 4, 8, 8, 8, 8, 12, 12, 12, 12, 16, 16, 16, 16/)

  integer, parameter :: &
    kfx = 300

  real(kind=8), parameter :: &
    zer = 0.0d+00, &
    pt5 = 0.5d+00, &
    one = 1.0d+00, &
    two = 2.0d+00, &
    ten = 10.0d+00, &
    pi = 4*atan(1.0_dp), &
    sqrtpi = sqrt(pi), &
    halfsqrtpi = 0.5d0*sqrtpi

  real(kind=8), parameter :: &
    tol = 1.0d-12

  real(kind=qp), parameter :: &
    oneq = 1.0e+00, &
    tolq = 1.0e-17

  integer :: max_m = 0

  real(kind=8) :: &
    fgrid(0:ntx_al, 0:npf, 0:ngrd), & !< interpolation tables for gamma
    xgrid(0:ntx_al, 0:npx), &        !< interpolation table for exp(-x)
    rfinc(0:ngrd), &             !< reciprocal interval widths for gamma
    rmr(mxqt), &                 !< 1/(2m-1) factors for fm(t) recursion
    tlgm(0:mxqt), &              !< (2m-1)!! factors for large-t formula
    tmax, &                      !< maximum value to interpolate for
    rxinc                        !< reciprocal interval width for exp(-x)

  integer :: &
      nord                         !< order of interpolation

!< Chebyshev polynomial coefficients (1st kind) up to 12th order (!)
!< Only nonzero coefficients, starting from the lowest power:
! GDF 12/09/02: do not know if 12th order works, tried 8th order.
  integer, parameter :: icheb(49) = &
      [ [ 1 ] &
      , [ 1 ] &
      , [ -1, 2 ] &
      , [ -3, 4 ] &
      , [  1, -8, 8 ] &
      , [  5, -20, 16 ] &
      , [ -1, 18, -48, 32 ] &
      , [ -7, 56, -112, 64 ] &
      , [  1, -32, 160, -256, 128 ] &
      , [  9, -120, 432, -576, 256 ] &
      , [ -1, 50, -400, 1120, -1280, 512 ] &
      , [ -11, 220, -1232, 2816, -2816, 1024 ] &
      , [  1, -72, 840, -3584, 6912, -6144, 2048 ] &
      ]
    !<  binomial coefficients - arranged to match Chebyshev terms
      integer, parameter :: ibino(252) = &
      [1,1,1,1,1,2,1,1,1,1,3,3,1,1,1,2,1,1,4,6,4,1            &
     ,          1,1,1,3,3,1,1,5,10,10,5,1,1,1,2,1,1,4,6,4,1            &
     ,          1,6,15,20,15,6,1,1,1,1,3,3,1,1,5,10,10,5,1,1           &
     ,          7,21,35,35,21,7,1,1,1,2,1,1,4,6,4,1,1,6,15,20          &
     ,          15,6,1,1,8,28,56,70,56,28,8,1,1,1,1,3,3,1,1,5          &
     ,          10,10,5,1,1,7,21,35,35,21,7,1,1,9,36,84,126,126        &
     ,          84,36,9,1,1,1,2,1,1,4,6,4,1,1,6,15,20,15,6,1,1,8       &
     ,          28,56,70,56,28,8,1,1,10,45,120,210,252,210,120,45      &
     ,          10,1,1,1,1,3,3,1,1,5,10,10,5,1,1,7,21,35,35,21,7,1     &
     ,          1,9,36,84,126,126,84,36,9,1,1,11,55,165,330,462,462    &
     ,          330,165,55,11,1,1,1,2,1,1,4,6,4,1,1,6,15,20,15,6,1,1   &
     ,          8,28,56,70,56,28,8,1,1,10,45,120,210,252,210,120,45    &
     ,          10,1,1,12,66,220,495,792,924,792,495,220,66,12,1]
contains

  subroutine boyspre(mmax)
!>    @brief   Boys function calculation pre-initialization
!>    @details Set up integral computation data for largest case.
!>             This routine can only be called when the highest
!>             angular momentum (l) of the basis has been established.
!>    @author  Graham Fletcher, 2004

    implicit none

    integer, intent(in) :: &
      mmax

    integer :: &
      ntms, ntol

    real(kind=8) :: &
      xmax

    integer :: &
      i, j, m

    real(kind=8) :: &
      fii

    ! generate interpolation grids
    ! this adds a couple of seconds to the total execution time,
    ! and requires the presence of quadruple precision.

    if (mmax <= max_m) return

    max_m = mmax

    ntms = ntx                 ! 4th order polynomial
    ntol = 12                ! 10**-12 accuracy
    xmax = 25.0d+00
    call cheby(mmax, ntms, ntol, xmax)

    ! reciprocal factors used in the downwards recursion for fm(t)

    do m = 1, mxqt
      rmr(m) = one/(2*m-1)
    end do

    ! (2m-1)!! factors for large-t (in fm(t))

    do i = 0, mmax
      fii = one
      do j = 1, 2*i-1, 2
        fii = j*fii
      end do
      tlgm(i) = fii
    end do

  end subroutine boyspre

  subroutine cheby(mmax, ntms, ntol, xmax)
!>    @brief   Incomplete gamma function interpolation setup
!>    @details This routine generates interpolation grids for the
!>             incomplete gamma function and exp(-x) based upon a
!>             Chebyshev polynomial fit.  In this approach the grid
!>             spacing is a function of the requested accuracy, so
!>             higher order (see ntms, below) corresponds to fewer
!>             grid points but more interpolation flops in the integral
!>             code, and vice versa for lower order. Gamma functions
!>             are computed using the Taylor series formula for small
!>             argument. Quad precision seems essential for obtaining
!>             accurate grid data near to the limit (25), especially
!>             for high order (see mmax, below).  A table of reciprocals
!>             is pre-computed for speed. for mmax higher than 16,
!>             kfx (below) may need to be increased.
!>    @author  Graham Fletcher, 2004

    implicit none

    integer, intent(in) :: &
      mmax, & !< highest order of gamma function
      ntms, & !< order of polynomial (4th is usually sufficient)
      ntol    !< accuracy of interpolation

    real(kind=8), intent(in) :: &
      xmax    !< end of interpolation range (25.0)

    real(kind=8) :: &
      xpw(0:ntx)

    real(kind=16) :: &
      rmk(kfx), qi

    real(kind=8) :: &
      tol, epn, &
      del, finc, xpt, xd, xpk, &
      cheb, awt, xinc

    integer :: &
      i, j, igrd, ipt, ktm, &
      nt2, ntf, ngrids, mft, npts, &
      ic, ib, ncheb, in, ii
!
    nord = ntms
    tmax = xmax
!
    tol = ten**(-ntol)
    nt2 = 2**ntms
    ntf = 1
    do i = 1, ntms+1
      ntf = ntf*i
    end do
    epn = one/(ntms+1)

    ! pre-compute reciprocals in quad

    do i = 1, kfx
      qi = i
      rmk(i) = 1.0_qp/qi
    end do

    ! Loop over selected gamma function orders, m
    ! grids are computed for m = 0,1,2,3,4, 8, 12, 16, ...
    ! To save storage, the grids for 8,12,16 are stored at fgrid(,,5),
    ! fgrid(,,6),fgrid(,,7) respectively, and you must recur downward
    ! from 8,12,16 wasting 8,7,6 values, if for example you want m=5.

    ngrids = 4
    if (mmax > 4) ngrids = ngrids+(mmax-1)/4
    do igrd = 0, ngrids
      mft = igrd
      if (igrd > 4) mft = (igrd-3)*4

      !  formula for optimum interval width

      del = nt2*ntf*(2*(mft+ntms+1)+1)
      del = del**epn
      del = del*(tol**epn)
      finc = del*two
      rfinc(igrd) = one/finc
      npts = nint(xmax/finc)

      ! loop over grid points

      xpt = zer
      do ipt = 0, npts

        ! generate powers of normalized grid point

        xd = -xpt/del
        xpk = one
        do j = 0, ntms
          xpw(j) = xpk
          xpk = xpk*xd
          fgrid(j, ipt, igrd) = zer
        end do

        ! generate term coefficients

        ic = 0
        ib = 0
        do ktm = 0, ntms
          call chebyg(del, xpt, ktm, mft, rmk, awt)
          ncheb = (ktm+2)/2    ! truncate
          in = mod(ktm, 2)
          do i = 1, ncheb
            ic = ic+1
            cheb = icheb(ic)*awt
            do j = 0, in
              ib = ib+1
              ii = in-j
              fgrid(ii, ipt, igrd) = fgrid(ii, ipt, igrd)+ &
                                            ibino(ib)*cheb*xpw(j)*(2**ii)
            end do
            in = in+2
          end do
        end do
        xpt = xpt+finc
      end do
    end do

    ! interpolation grid for exp(-x)
    ! formula for optimum interval width

    del = nt2*ntf
    del = del**epn
    del = del*(tol**epn)
    xinc = del*two
    rxinc = one/xinc
    npts = nint(xmax/xinc)

    ! loop over grid points

    xpt = zer
    do ipt = 0, npts

      ! generate powers of normalized grid point
      xd = -xpt/del
      xpk = one
      do j = 0, ntms
        xpw(j) = xpk
        xpk = xpk*xd
        xgrid(j, ipt) = zer
      end do

      !  generate term coefficients

      ic = 0
      ib = 0
      do ktm = 0, ntms
        call chebyx(del, xpt, ktm, awt)
        ncheb = (ktm+2)/2    ! truncate
        in = mod(ktm, 2)
        do i = 1, ncheb
          ic = ic+1
          cheb = icheb(ic)*awt
          do j = 0, in
            ib = ib+1
            ii = in-j
            xgrid(ii, ipt) = xgrid(ii, ipt)+ &
                                    ibino(ib)*cheb*xpw(j)*(2**ii)
          end do
          in = in+2
        end do
      end do
      xpt = xpt+xinc
    end do

  end subroutine cheby

  subroutine chebyg(del, xpt, ktm, mft, rmk, awt)
!>    @brief   Chebyshev weight for gamma function
!>    @details Chebyshev weight for gamma function
!>    @author  Graham Fletcher, 2004

    implicit none

    real(kind=8), intent(in) :: &
      del, &    !< half interval width
      xpt       !< grid point

    integer, intent(in) :: &
      ktm, &    !< polynomial term
      mft       !< order of gamma function

    real(kind=16), intent(in) :: &
      rmk(kfx)  !< table of reciprocals
    !< quad not essential here but hey

    real(kind=8), intent(out) :: &
      awt       !< term weight

    logical :: &
      converged

    real(kind=8) :: &
      dd2, dsq, am, fk, fm, ap, ftm

    integer :: &
      m, mmk

    dd2 = del*pt5
    dsq = dd2*dd2
    awt = zer
    am = one
    fk = one
    do m = 1, ktm
      fk = real(fk*rmk(m), dp)
    end do

    converged = .false.
    m = 0
    fm = one
    do while (.not. converged)
      ap = awt
      mmk = ktm+mft+2*m
      call ftmval(xpt, mmk, rmk, ftm)
      awt = awt+am*ftm*fm*fk
      converged = abs(ap-awt) < tol
      m = m+1
      fm = real(fm*rmk(m), dp)
      fk = real(fk*rmk(ktm+m), dp)
      am = am*dsq
    end do
    awt = awt*(dd2**ktm)
    if (ktm > 0) awt = awt*two
    if (mod(ktm, 2) /= 0) awt = -awt

  end subroutine chebyg

  subroutine ftmval(tt8, mft, rmk, ftm8)
!>    @brief   evaluate Fm(t) table
!>    @details Formula for computing the gamma function using Taylor
!>             expansion for tt<25. care taken to do this in a
!>             numerically stable way.  Note that quadruple precision
!>             is used, although double precision arguments are passed.
!>             quadruple precision is critical!
!>             Because Q.P. is frequently unavailable, this routine
!>             is not called, instead the data it generates is read
!>             from a disk file prepared by a machine that has Q.P.
!>    @author  Graham Fletcher, 2004

    implicit none

    real(kind=8), intent(in) :: &
      tt8         !< gamma function argument

    real(kind=16), intent(in) :: &
      rmk(kfx)    !< table of reciprocals

    integer, intent(in) :: &
      mft         !< order of gamma function

    real(kind=8), intent(out) :: &
      ftm8        !< value of gamma function

    logical :: &
      converged

    real(kind=16) :: &
      tt, ftm, xk, a, fp

    integer :: &
      l, m

    tt = tt8

    converged = .false.
    l = 1
    m = 2*mft+1
    ftm = rmk(m)
    xk = oneq
    do while (.not. converged)
      m = m+2
      a = -tt*rmk(m)
      fp = ftm
      ftm = ftm+a*xk
      converged = abs(fp-ftm) < tolq
      l = l+1
      xk = -xk*tt*rmk(l)
    end do

    ftm8 = real(ftm, dp)

  end subroutine ftmval

  subroutine chebyx(del, xpt, ktm, awt)
!>    @brief   Chebyshev weight for exp(-x)
!>    @details Chebyshev weight for exp(-x)
!>             Note that apart from the last line the weights
!>             are independent of the value of xpt
!>    @author  Graham Fletcher, 2004

    implicit double precision(a-h, o-z)

    integer, intent(in) :: &
      ktm    !< polynomial term

    real(kind=8), intent(in) :: &
      del, & !< half interval width
      xpt    !< grid point

    real(kind=8), intent(out) :: &
      awt    !< term weight

    logical :: &
      converged

    real(kind=8) :: &
      dd2, dsq, am, fk, fm, ap

    integer :: &
      m

    dd2 = del*pt5
    dsq = dd2*dd2
    awt = zer
    am = one
    fk = one
    do m = 1, ktm
      fk = fk*m
    end do

    converged = .false.
    m = 0
    fm = one
    do while (.not. converged)
      ap = awt
      awt = awt+am/(fm*fk)
      converged = abs(ap-awt) < tol
      m = m+1
      fm = fm*m
      fk = fk*(ktm+m)
      am = am*dsq
    end do
    awt = awt*(dd2**ktm)
    if (ktm > 0) awt = awt*two
    if (mod(ktm, 2) /= 0) awt = -awt
    awt = awt*exp(-xpt)

  end subroutine chebyx
end module

program gen_boys
  use boys_prepare
  integer, parameter :: i = 16
  character(*), parameter :: int_name = 'integer, public, parameter'
  character(*), parameter :: real_name = 'real(kind=8), public, parameter'
  character(*), parameter :: int_fmt = '(2XA/ *( 4x, 20(I0,:,", "), "&", /) )'
  character(*), parameter :: real_fmt = '(2XA/ *( 4x, 3(D23.16,:,", "), "&", /) )'

  call boyspre(i)
  write(*,'(G0)') 'module boys_lut'
  write(*,'(2XA)') 'implicit none'
  write(*,'(2XA)') 'private'

  write(*,'(2XA,I0)') int_name//' :: ntx  = ', ntx
  write(*,'(2XA,I0)') int_name//' :: npf  = ', npf
  write(*,'(2XA,I0)') int_name//' :: ngrd = ', ngrd
  write(*,'(2XA,I0)') int_name//' :: npx  = ', npx
  write(*,'(2XA,I0)') int_name//' :: mxqt = ', mxqt
  write(*,'(2XA,I0)') int_name//' :: ntx_al = ', ntx_al

  write(*,advance='no',fmt=int_fmt) int_name//' :: igrid(0:16) = [ &', igrid
  write(*,*) ']'
  write(*,advance='no',fmt=int_fmt) int_name//' :: irgrd(0:16) = [ &', igrid
  write(*,*) ']'

  write(*,advance='no',fmt=real_fmt) real_name//' :: fgrid(0:ntx_al,0:npf,0:ngrd) = reshape([ &', fgrid
  write(*,*) '], shape=[ntx_al+1, npf+1, ngrd+1])'

  write(*,advance='no',fmt=real_fmt) real_name//' :: xgrid(0:ntx_al,0:npx) = reshape([ &', xgrid
  write(*,*) '], shape=[ntx_al+1, npx+1])'

  write(*,advance='no',fmt=real_fmt) real_name//' :: rfinc(0:ngrd) = [ &', rfinc
  write(*,*) ']'

  write(*,advance='no',fmt=real_fmt) real_name//' :: rmr(mxqt) = [ &', rmr
  write(*,*) ']'

  write(*,advance='no',fmt=real_fmt) real_name//' :: tlgm(0:mxqt) = [ &', tlgm
  write(*,*) ']'

  write(*,'(2XA,D23.16)') real_name//' :: tmax = ', tmax

  write(*,'(2XA,D23.16)') real_name//' :: rxinc = ', rxinc

  write(*,'(2XA,I0)') int_name//' :: nord = ', nord

  write(*,'(G0)') 'end module boys_lut'
end program gen_boys
