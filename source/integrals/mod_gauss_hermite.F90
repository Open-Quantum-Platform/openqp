!#define DEBUG
! Protecting macro for compilers which does not support OpenMP 4.0
#define OMPSIMD (_OPENMP >= 201307)

!> @brief Gauss-Hermite quadrature used in one-electron integral code
!> @author   Vladimir Mironov
!
!     REVISION HISTORY:
!> @date _Sep, 2018_ Initial release
!
 module mod_gauss_hermite

    use precision, only: dp

!    use constants, only: h => hermit, w => hermitw
    implicit none

    private
    public doQuadGaussHermite
    public derQuadGaussHermite
    public mulQuadGaussHermite
!  Roots and weights for Gauss-Hermite quadrature
!  The values used below were obtained by running the routine given
!  in the "numerical recipes" book in quadruple precision.
    real(kind=dp), parameter :: h2d(10,10)=reshape([&
       0.0000000000000000000D+00,  0.0000000000000000000D+00,  0.0000000000000000000D+00,  0.0000000000000000000D+00,  &
       0.0000000000000000000D+00,  0.0000000000000000000D+00,  0.0000000000000000000D+00,  0.0000000000000000000D+00,  &
       0.0000000000000000000D+00,  0.0000000000000000000D+00, -7.0710678118654752440D-01,  7.0710678118654752440D-01,  &
       0.0000000000000000000D+00,  0.0000000000000000000D+00, 0.0000000000000000000D+00,  0.0000000000000000000D+00,   &
       0.0000000000000000000D+00,  0.0000000000000000000D+00, 0.0000000000000000000D+00,  0.0000000000000000000D+00,   &
       -1.2247448713915890491D+00,  0.0000000000000000000D+00,  1.2247448713915890491D+00,  0.0000000000000000000D+00, &
       0.0000000000000000000D+00,  0.0000000000000000000D+00,  0.0000000000000000000D+00,  0.0000000000000000000D+00,  &
       0.0000000000000000000D+00,  0.0000000000000000000D+00, -1.6506801238857845559D+00, -5.2464762327529031788D-01,  &
       5.2464762327529031788D-01,  1.6506801238857845559D+00, 0.0000000000000000000D+00,  0.0000000000000000000D+00,   &
       0.0000000000000000000D+00,  0.0000000000000000000D+00, 0.0000000000000000000D+00,  0.0000000000000000000D+00,   &
       -2.0201828704560856329D+00, -9.5857246461381850711D-01,  0.0000000000000000000D+00,  9.5857246461381850711D-01, &
       2.0201828704560856329D+00,  0.0000000000000000000D+00,  0.0000000000000000000D+00,  0.0000000000000000000D+00,  &
       0.0000000000000000000D+00,  0.0000000000000000000D+00, -2.3506049736744922228D+00, -1.3358490740136969497D+00,  &
       -4.3607741192761650868D-01,  4.3607741192761650868D-01, 1.3358490740136969497D+00,  2.3506049736744922228D+00,  &
       0.0000000000000000000D+00,  0.0000000000000000000D+00, 0.0000000000000000000D+00,  0.0000000000000000000D+00,   &
       -2.6519613568352334925D+00, -1.6735516287674714450D+00, -8.1628788285896466304D-01,  0.0000000000000000000D+00, &
       8.1628788285896466304D-01,  1.6735516287674714450D+00,  2.6519613568352334925D+00,  0.0000000000000000000D+00,  &
       0.0000000000000000000D+00,  0.0000000000000000000D+00, -2.9306374202572440192D+00, -1.9816567566958429259D+00,  &
       -1.1571937124467801947D+00, -3.8118699020732211685D-01, 3.8118699020732211685D-01,  1.1571937124467801947D+00,  &
       1.9816567566958429259D+00,  2.9306374202572440192D+00, 0.0000000000000000000D+00,  0.0000000000000000000D+00,   &
       -3.1909932017815276072D+00, -2.2665805845318431118D+00, -1.4685532892166679317D+00, -7.2355101875283757332D-01, &
       0.0000000000000000000D+00,  7.2355101875283757332D-01,  1.4685532892166679317D+00,  2.2665805845318431118D+00,  &
       3.1909932017815276072D+00,  0.0000000000000000000D+00, -3.4361591188377376033D+00, -2.5327316742327897964D+00,  &
       -1.7566836492998817735D+00, -1.0366108297895136542D+00, -3.4290132722370460879D-01,  3.4290132722370460879D-01, &
       1.0366108297895136542D+00,  1.7566836492998817735D+00, 2.5327316742327897964D+00,  3.4361591188377376033D+00  &
       ], shape=shape(h2d))
    real(kind=dp), parameter :: w2d(10,10)=reshape([&
       1.7724538509055160273D+00,  0.0000000000000000000D+00,  0.0000000000000000000D+00,  0.0000000000000000000D+00, &
       0.0000000000000000000D+00,  0.0000000000000000000D+00,  0.0000000000000000000D+00,  0.0000000000000000000D+00, &
       0.0000000000000000000D+00,  0.0000000000000000000D+00, 8.8622692545275801365D-01,  8.8622692545275801365D-01, &
       0.0000000000000000000D+00,  0.0000000000000000000D+00, 0.0000000000000000000D+00,  0.0000000000000000000D+00, &
       0.0000000000000000000D+00,  0.0000000000000000000D+00, 0.0000000000000000000D+00,  0.0000000000000000000D+00, &
       2.9540897515091933788D-01,  1.1816359006036773515D+00,  2.9540897515091933788D-01,  0.0000000000000000000D+00, &
       0.0000000000000000000D+00,  0.0000000000000000000D+00,  0.0000000000000000000D+00,  0.0000000000000000000D+00, &
       0.0000000000000000000D+00,  0.0000000000000000000D+00, 8.1312835447245177143D-02,  8.0491409000551283651D-01, &
       8.0491409000551283651D-01,  8.1312835447245177143D-02, 0.0000000000000000000D+00,  0.0000000000000000000D+00, &
       0.0000000000000000000D+00,  0.0000000000000000000D+00, 0.0000000000000000000D+00,  0.0000000000000000000D+00, &
       1.9953242059045913208D-02,  3.9361932315224115983D-01,  9.4530872048294188123D-01,  3.9361932315224115983D-01, &
       1.9953242059045913208D-02,  0.0000000000000000000D+00,  0.0000000000000000000D+00,  0.0000000000000000000D+00, &
       0.0000000000000000000D+00,  0.0000000000000000000D+00, 4.5300099055088456409D-03,  1.5706732032285664392D-01, &
       7.2462959522439252409D-01,  7.2462959522439252409D-01, 1.5706732032285664392D-01,  4.5300099055088456409D-03, &
       0.0000000000000000000D+00,  0.0000000000000000000D+00, 0.0000000000000000000D+00,  0.0000000000000000000D+00, &
       9.7178124509951915415D-04,  5.4515582819127030592D-02,  4.2560725261012780052D-01,  8.1026461755680732676D-01, &
       4.2560725261012780052D-01,  5.4515582819127030592D-02,  9.7178124509951915415D-04,  0.0000000000000000000D+00, &
       0.0000000000000000000D+00,  0.0000000000000000000D+00, 1.9960407221136761921D-04,  1.7077983007413475456D-02, &
       2.0780232581489187954D-01,  6.6114701255824129103D-01, 6.6114701255824129103D-01,  2.0780232581489187954D-01, &
       1.7077983007413475456D-02,  1.9960407221136761921D-04, 0.0000000000000000000D+00,  0.0000000000000000000D+00, &
       3.9606977263264381905D-05,  4.9436242755369472172D-03,  8.8474527394376573288D-02,  4.3265155900255575020D-01, &
       7.2023521560605095712D-01,  4.3265155900255575020D-01,  8.8474527394376573288D-02,  4.9436242755369472172D-03, &
       3.9606977263264381905D-05,  0.0000000000000000000D+00, 7.6404328552326206292D-06,  1.3436457467812326922D-03, &
       3.3874394455481063136D-02,  2.4013861108231468642D-01, 6.1086263373532579878D-01,  6.1086263373532579878D-01, &
       2.4013861108231468642D-01,  3.3874394455481063136D-02, 1.3436457467812326922D-03,  7.6404328552326206292D-06 &
       ], shape=shape(w2d))

contains

!--------------------------------------------------------------------------------

!> @brief Gauss-Hermite quadrature using minimum point formula
!> @details Compute:
!>  xint = sum( w(1:npts,npts) * (h(1:npts,npts)*t+dxi)**(ni-1) * (h(1:npts,npts)*t+dxj)**(nj-1) )
!>  yint = sum( w(1:npts,npts) * (h(1:npts,npts)*t+dyi)**(ni-1) * (h(1:npts,npts)*t+dyj)**(nj-1) )
!>  zint = sum( w(1:npts,npts) * (h(1:npts,npts)*t+dzi)**(ni-1) * (h(1:npts,npts)*t+dzj)**(nj-1) )
!> @note Use of I functions will let NI run up to 7 (S=1, P=2, ..., I=7)
!>       Use of I functions will let NJ run up to 9 (for K.E. ints)
!>       Use of I functions requires NPTS=8 to do kinetic energy integrals.
!> @param[out]      xint        x-component of the integral
!> @param[out]      yint        y-component of the integral
!> @param[out]      zint        z-component of the integral
!> @param[out]      t           inverse square root of total exponent
!> @param[out]      x0          `x`-coord. of the center of primitive pair
!> @param[out]      y0          `y`-coord. of the center of primitive pair
!> @param[out]      z0          `z`-coord. of the center of primitive pair
!> @param[out]      xi          `x`-coord. of the center of 1st primitive
!> @param[out]      yi          `y`-coord. of the center of 1st primitive
!> @param[out]      zi          `z`-coord. of the center of 1st primitive
!> @param[out]      xj          `x`-coord. of the center of 2nd primitive
!> @param[out]      yj          `y`-coord. of the center of 2nd primitive
!> @param[out]      zj          `z`-coord. of the center of 2nd primitive
!> @param[out]      ni          current angular momentum on center i
!> @param[out]      nj          current angular momentum on center j
!
!> @note based on STVINT from int1.src
!
!> @author   Vladimir Mironov
!
!     REVISION HISTORY:
!> @date _Sep, 2018_ Initial release
!
 subroutine doQuadGaussHermite(tint, t, rij, ri, rj, ni, nj)
!dir$ attributes forceinline :: doQuadGaussHermite

 real(kind=dp), intent(out) :: tint(3)
 real(kind=dp), intent(in)  :: t
 real(kind=dp), intent(in)  :: rij(3), ri(3), rj(3)
 integer, intent(in) :: ni, nj

 real(kind=dp) :: dqi(3), dqj(3), p(3)
 integer :: i, j, npts

#ifdef DEBUG
    if ( ni>7 .or. nj>9 .or. ni<1 .or. nj<1 ) then
        write (*,'(" stvint: exceeded limitations, with ni,nj=",2i5)') ni, nj
        call abrt
        return
    end if
#endif

    npts = (ni+nj-2)/2 + 1

    dqi = rij - ri

    dqj = rij - rj

    tint = 0.0

#if OMPSIMD
!!$omp simd reduction(+:tint) &
!!$omp   private(p)
#endif
    do i = 1, npts

        p = w2d(i,npts)

        do j = 2, ni
            p = p*(h2d(i,npts)*t + dqi)
        end do

        do j = 2, nj
            p = p*(h2d(i,npts)*t + dqj)
        end do

        tint = tint + p

    end do

 end subroutine

!--------------------------------------------------------------------------------

!> @brief Gauss-Hermite quadrature using minimum point formula. Coulomb 1e derivatives variant.
!> @details Compute:
!> `xint = sum( (h(1:npts,npts)*t+dxc) * w(1:npts,npts) * (h(1:npts,npts)*t+dxi)**(ni-1) * (h(1:npts,npts)*t+dxj)**(nj-1) )`
!> `yint = sum( (h(1:npts,npts)*t+dyc) * w(1:npts,npts) * (h(1:npts,npts)*t+dyi)**(ni-1) * (h(1:npts,npts)*t+dyj)**(nj-1) )`
!> `zint = sum( (h(1:npts,npts)*t+dzc) * w(1:npts,npts) * (h(1:npts,npts)*t+dzi)**(ni-1) * (h(1:npts,npts)*t+dzj)**(nj-1) )`
!
!> @param[out]      xint        x-component of the integral
!> @param[out]      yint        y-component of the integral
!> @param[out]      zint        z-component of the integral
!> @param[out]      t           inverse square root of total exponent
!> @param[out]      x0          `x`-coord. of the center of primitive pair
!> @param[out]      y0          `y`-coord. of the center of primitive pair
!> @param[out]      z0          `z`-coord. of the center of primitive pair
!> @param[out]      xi          `x`-coord. of the center of 1st primitive
!> @param[out]      yi          `y`-coord. of the center of 1st primitive
!> @param[out]      zi          `z`-coord. of the center of 1st primitive
!> @param[out]      xj          `x`-coord. of the center of 2nd primitive
!> @param[out]      yj          `y`-coord. of the center of 2nd primitive
!> @param[out]      zj          `z`-coord. of the center of 2nd primitive
!> @param[out]      cx          `x`-coord. of the charged particle
!> @param[out]      cy          `y`-coord. of the charged particle
!> @param[out]      cz          `z`-coord. of the charged particle
!> @param[out]      ni          current angular momentum on center i
!> @param[out]      nj          current angular momentum on center j
!
!> @note based on DVINT from grd1.src
!> @note currently not used, left for debugging purposes
!
!> @author   Vladimir Mironov
!
!     REVISION HISTORY:
!> @date _Sep, 2018_ Initial release
!
 subroutine derQuadGaussHermite(xint, yint, zint, t, x0, y0, z0, xi, yi, zi, xj, yj, zj, cx, cy, cz, ni, nj)
!dir$ attributes forceinline :: derQuadGaussHermite

 real(kind=dp), intent(out) :: xint, yint, zint
 real(kind=dp), intent(in)  :: t, x0, y0, z0, xi, yi, zi, xj, yj, zj, cx, cy, cz
 integer, intent(in) :: ni, nj

 real(kind=dp) :: dxi, dyi, dzi, dxj, dyj, dzj, dxc, dyc, dzc, px, py, pz
 integer :: i, j, npts

    npts = (ni+nj-2+1)/2 + 1

    dxi = x0 - xi
    dyi = y0 - yi
    dzi = z0 - zi

    dxj = x0 - xj
    dyj = y0 - yj
    dzj = z0 - zj

    dxc = x0 - cx
    dyc = y0 - cy
    dzc = z0 - cz

    xint = 0.0
    yint = 0.0
    zint = 0.0

#if OMPSIMD
!$omp simd reduction(+:xint,yint,zint) &
!$omp   private(px,py,pz)
#endif
    do i = 1, npts

        px = h2d(i,npts)*t + dxc
        py = h2d(i,npts)*t + dyc
        pz = h2d(i,npts)*t + dzc

        do j = 2, ni
            px = px*(h2d(i,npts)*t + dxi)
            py = py*(h2d(i,npts)*t + dyi)
            pz = pz*(h2d(i,npts)*t + dzi)
        end do

        do j = 2, nj
            px = px*(h2d(i,npts)*t + dxj)
            py = py*(h2d(i,npts)*t + dyj)
            pz = pz*(h2d(i,npts)*t + dzj)
        end do

        xint = xint + w2d(i,npts)*px
        yint = yint + w2d(i,npts)*py
        zint = zint + w2d(i,npts)*pz
    end do

 end subroutine

!--------------------------------------------------------------------------------

!> @brief Gauss-Hermite quadrature using minimum point formula
!> @details Compute:
!>  xint = sum( w(1:npts,npts) * (h(1:npts,npts)*t+dxi)**(ni-1) * (h(1:npts,npts)*t+dxj)**(nj-1) )
!>  yint = sum( w(1:npts,npts) * (h(1:npts,npts)*t+dyi)**(ni-1) * (h(1:npts,npts)*t+dyj)**(nj-1) )
!>  zint = sum( w(1:npts,npts) * (h(1:npts,npts)*t+dzi)**(ni-1) * (h(1:npts,npts)*t+dzj)**(nj-1) )
!> @note Use of I functions will let NI run up to 7 (S=1, P=2, ..., I=7)
!>       Use of I functions will let NJ run up to 9 (for K.E. ints)
!>       Use of I functions requires NPTS=8 to do kinetic energy integrals.
!> @param[out]      xint        x-component of the integral
!> @param[out]      yint        y-component of the integral
!> @param[out]      zint        z-component of the integral
!> @param[out]      t           inverse square root of total exponent
!> @param[out]      x0          `x`-coord. of the center of primitive pair
!> @param[out]      y0          `y`-coord. of the center of primitive pair
!> @param[out]      z0          `z`-coord. of the center of primitive pair
!> @param[out]      xi          `x`-coord. of the center of 1st primitive
!> @param[out]      yi          `y`-coord. of the center of 1st primitive
!> @param[out]      zi          `z`-coord. of the center of 1st primitive
!> @param[out]      xj          `x`-coord. of the center of 2nd primitive
!> @param[out]      yj          `y`-coord. of the center of 2nd primitive
!> @param[out]      zj          `z`-coord. of the center of 2nd primitive
!> @param[out]      ni          current angular momentum on center i
!> @param[out]      nj          current angular momentum on center j
!
!> @note based on STVINT from int1.src
!
!> @author   Vladimir Mironov
!
!     REVISION HISTORY:
!> @date _Sep, 2018_ Initial release
!
 subroutine mulQuadGaussHermite(tint, t, rij, ri, rj, r, ni, nj, m)
!dir$ attributes forceinline :: doQuadGaussHermite

 real(kind=dp), intent(out) :: tint(3,0:*)
 real(kind=dp), intent(in)  :: t
 real(kind=dp), intent(in)  :: rij(3), ri(3), rj(3), r(3)
 integer, intent(in) :: ni, nj, m

 real(kind=dp) :: dqi(3), dqj(3), dqc(3), p(3)
 integer :: i, j, npts

#ifdef DEBUG
    if ( ni>7 .or. nj>9 .or. ni<1 .or. nj<1 ) then
        write (*,'(" stvint: exceeded limitations, with ni,nj=",2i5)') ni, nj
        call abrt
        return
    end if
#endif

    npts = (ni+nj+m-2)/2 + 1

    dqi = rij - ri

    dqj = rij - rj

    dqc = rij - r

    tint(:,0:m) = 0.0

#if OMPSIMD
!!$omp simd reduction(+:tint) &
!!$omp   private(p)
#endif
    do i = 1, npts

        p = w2d(i,npts)

        do j = 2, ni
            p = p*(h2d(i,npts)*t + dqi)
        end do

        do j = 2, nj
            p = p*(h2d(i,npts)*t + dqj)
        end do

        do j = 0, m
            tint(:,j) = tint(:,j) + p
            p = p*(h2d(i,npts)*t + dqc)
        end do

    end do

 end subroutine

 end module
