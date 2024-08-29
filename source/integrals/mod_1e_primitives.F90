! Protecting macro for compilers which does not support OpenMP 4.0
!#define OMPSIMD (_OPENMP >= 201307)

!> @brief Helper functions and data blocks needed
!> to compute one-electron integrals and their derivatives
!
!> @author   Vladimir Mironov
!
!> @todo
!> - Unify interfaces
!> - Cleanup redundant subroutines
!
!     REVISION HISTORY:
!> @date _Sep, 2018_ Initial release
!
MODULE mod_1e_primitives
 USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: REAL64
 USE mod_gauss_hermite, ONLY: doQuadGaussHermite, mulQuadGaussHermite
 USE mod_shell_tools, ONLY: shell_t, shpair_t
 use rys, only: rys_root_t
 use xyz_order
 use constants, only: PI, CART_X, CART_Y, CART_Z, MAX_ANG => BAS_MXANG
 IMPLICIT NONE

 INTEGER :: iii
 !integer, parameter :: MAX_ANG = 6
 integer, parameter :: MAX_ANG_PAD = 7
 integer, parameter :: MAX_NROOTS = (2*MAX_ANG+1)/2+1
 integer, parameter, public :: MAX_EL_MOM = 3
 character, parameter :: MAX_EL_MOM_S = '3'

 REAL(REAL64), PARAMETER :: TWOPI  = pi * 2.0_real64


 PRIVATE

 PUBLIC comp_coulomb_int1_prim
 PUBLIC comp_kin_ovl_int1_prim
 PUBLIC comp_lz_int1_prim
 public comp_mult_int1_prim
 public comp_allmult_int1_prim
 PUBLIC comp_coulomb_dampch_int1_prim
 PUBLIC comp_ewaldlr_int1_prim
 PUBLIC comp_coulpot_prim
 PUBLIC comp_coulomb_der1
 PUBLIC comp_coulomb_helfeyder1
 PUBLIC comp_kinetic_der1
 PUBLIC comp_overlap_der1
 PUBLIC comp_ewaldlr_der1
 PUBLIC comp_ewaldlr_helfeyder1

 PUBLIC update_triang_matrix
 PUBLIC update_rectangular_matrix
 PUBLIC density_ordered
 PUBLIC density_unordered

CONTAINS


!--------------------------------------------------------------------------------
!       ONE-ELECTRON INTEGRALS CALCULATION (PRIMITIVE GAUSSIANS)
!--------------------------------------------------------------------------------

!> @brief Compute primitive block of overlap and kinetic energy 1e integrals
!> @param[in]       cp          shell pair data
!> @param[in]       id          current pair of primitives
!> @param[in]       dokinetic   if `.FALSE.` compute only overlap integrals
!> @param[inout]    sblk        block of 1e overlap integrals
!> @param[inout]    tblk        block of 1e kinetic energy integrals
!
!> @author   Vladimir Mironov
!
!     REVISION HISTORY:
!> @date _Sep, 2018_ Initial release
!
 SUBROUTINE comp_kin_ovl_int1_prim(cp, id, dokinetic, sblk, tblk)
!dir$ attributes inline :: comp_kin_ovl_int1_prim

    TYPE(shpair_t), INTENT(IN)  :: cp
    INTEGER, INTENT(IN) :: id
    LOGICAL, INTENT(IN)   :: dokinetic
    REAL(REAL64), CONTIGUOUS,   INTENT(INOUT) :: sblk(:), tblk(:)

    INTEGER :: i, j, nx, ny, nz, mx, my, mz, jmax, ij
    REAL(REAL64) :: ovl, kinx, kiny, kinz, kin
    real(real64) :: xyzovl(0:max_ang+2,0:max_ang,3)
    real(real64) :: xyzkin(0:max_ang_pad,0:max_ang,3)
!dir$ assume_aligned sblk : 64
!dir$ assume_aligned tblk : 64
!dir$ assume_aligned xyzkin : 64
!dir$ assume_aligned xyzovl : 64

    jmax = cp%jang
    IF (dokinetic) jmax = cp%jang+2

    ASSOCIATE (pp => cp%p(id))
    CALL overlap_xyz(cp%ri, cp%rj, pp%r, pp%aa1, cp%iang, jmax, xyzovl)

    IF (dokinetic) CALL kinetic_xyz_j(xyzkin, xyzovl, cp%iang, cp%jang, pp%aj)

    ij = 0
    jmax = cp%jnao
    DO i = 1, cp%inao
        nx = CART_X(i,cp%iang)
        ny = CART_Y(i,cp%iang)
        nz = CART_Z(i,cp%iang)
        IF (cp%iandj) jmax = i
        DO j = 1, jmax
            mx = CART_X(j,cp%jang)
            my = CART_Y(j,cp%jang)
            mz = CART_Z(j,cp%jang)
            ij = ij+1

            ovl = xyzovl(mx,nx,1)*xyzovl(my,ny,2)*xyzovl(mz,nz,3)
            sblk(ij) = sblk(ij) + pp%expfac*ovl

            IF (dokinetic) THEN
                kinx = xyzkin(mx,nx,1)*xyzovl(my,ny,2)*xyzovl(mz,nz,3)
                kiny = xyzovl(mx,nx,1)*xyzkin(my,ny,2)*xyzovl(mz,nz,3)
                kinz = xyzovl(mx,nx,1)*xyzovl(my,ny,2)*xyzkin(mz,nz,3)
                kin  = kinx + kiny + kinz
                tblk(ij) =  tblk(ij) + pp%expfac*kin
            END IF
        END DO

    END DO

    END ASSOCIATE

 END SUBROUTINE

!> @brief Compute primitive block of 1e Coulomb atraction integrals
!> @param[in]       cp          shell pair data
!> @param[in]       id          current pair of primitives
!> @param[in]       c           coordinates of the charged particle
!> @param[in]       znuc        particle charge
!> @param[inout]    vblk        block of 1e Coulomb integrals
!
!> @author   Vladimir Mironov
!
!     REVISION HISTORY:
!> @date _Sep, 2018_ Initial release
!
 SUBROUTINE comp_coulomb_int1_prim(cp, id, c, znuc, vblk)
!dir$ attributes inline :: comp_coulomb_int1_prim

    TYPE(shpair_t), INTENT(IN)    :: cp
    INTEGER, INTENT(IN)   :: id
    REAL(REAL64), INTENT(IN)   :: c(3), znuc
    REAL(REAL64), CONTIGUOUS,   INTENT(INOUT) :: vblk(:)

      REAL(REAL64) :: xx

    type(rys_root_t) :: ryscomp
    INTEGER :: i, j, ij, jmax, nx, ny, nz, mx, my, mz
    REAL(REAL64) :: dum, dij
    real(real64) :: xyzin(0:2*max_ang+1, 0:max_ang+1,3,max_nroots)
!dir$ assume_aligned xyzin : 64
!dir$ assume_aligned vblk : 64

    ASSOCIATE (pp => cp%p(id), &
               iang => cp%iang, jang => cp%jang, &
               inao => cp%inao, jnao => cp%jnao)

    xx = pp%aa* sum((pp%r-c)**2)
    ryscomp%nroots = cp%nroots
    ryscomp%x = xx
    CALL QGaussRys(ryscomp, cp, id, c, znuc, xyzin)

    dij = pp%expfac*TWOPI*pp%aa1
    ij = 0
    jmax = jnao
    DO i = 1, inao
        nx = CART_X(i,iang)
        ny = CART_Y(i,iang)
        nz = CART_Z(i,iang)
        IF (cp%iandj) jmax = i
        DO j = 1, jmax
            mx = CART_X(j,jang)
            my = CART_Y(j,jang)
            mz = CART_Z(j,jang)
            ij = ij+1

            dum = dij * sum(  xyzin(mx,nx,1,1:cp%nroots) &
                            * xyzin(my,ny,2,1:cp%nroots) &
                            * xyzin(mz,nz,3,1:cp%nroots) )

            vblk(ij) = vblk(ij) + dum
        END DO
    END DO
    END ASSOCIATE

END SUBROUTINE

!> @brief Compute primitive block of 1e Coulomb atraction integrals for
!>  Ewald summation, long-range part
!> @details 1e integrals using modified Coulomb potential:
!>  \f$ \frac{Erf(\omega^{1/2}|r-r_C|)}{|r-r_C|} \f$
!> @param[in]       cp          shell pair data
!> @param[in]       id          current pair of primitives
!> @param[in]       c           coordinates of the charged particle
!> @param[in]       znuc        particle charge
!> @param[in]       omega       Ewald splitting parameter
!> @param[inout]    vblk        block of 1e Coulomb integrals
!
!> @author   Vladimir Mironov
!
!     REVISION HISTORY:
!> @date _Sep, 2018_ Initial release
!
 SUBROUTINE comp_ewaldlr_int1_prim(cp, id, c, znuc, omega, vblk)
!dir$ attributes inline :: comp_ewaldlr_int1_prim

    TYPE(shpair_t), INTENT(IN)    :: cp
    INTEGER, INTENT(IN)   :: id
    REAL(REAL64), INTENT(IN)   :: c(3), znuc
    REAL(REAL64), INTENT(IN)   :: omega
    REAL(REAL64), CONTIGUOUS,   INTENT(INOUT) :: vblk(:)

    REAL(REAL64) :: xx
    type(rys_root_t) :: ryscomp

    INTEGER :: i, j, ij, jmax, nx, ny, nz, mx, my, mz
    REAL(REAL64) :: dum, dij, xfac
    real(real64) :: xyzin(0:2*max_ang+1, 0:max_ang+1,3,max_nroots)
!dir$ assume_aligned xyzin : 64
!dir$ assume_aligned vblk : 64

    ASSOCIATE (pp => cp%p(id), &
               iang => cp%iang, jang => cp%jang, &
               inao => cp%inao, jnao => cp%jnao)

    IF (omega<=0.0d0) RETURN

    xfac = omega*omega/(pp%aa+omega*omega)
    xx = pp%aa * sum((pp%r-c)**2) * xfac
    ryscomp%nroots = cp%nroots
    ryscomp%x = xx
    CALL QGaussRysEw(ryscomp, cp, id, c, znuc, xfac, xyzin)

    dij = pp%expfac*TWOPI*pp%aa1 * sqrt(xfac)
    ij = 0
    jmax = jnao
    DO i = 1, inao
        nx = CART_X(i,iang)
        ny = CART_Y(i,iang)
        nz = CART_Z(i,iang)
        IF (cp%iandj) jmax = i
        DO j = 1, jmax
            mx = CART_X(j,jang)
            my = CART_Y(j,jang)
            mz = CART_Z(j,jang)
            ij = ij+1

            dum = dij * sum(  xyzin(mx,nx,1,1:cp%nroots) &
                            * xyzin(my,ny,2,1:cp%nroots) &
                            * xyzin(mz,nz,3,1:cp%nroots) )

            vblk(ij) = vblk(ij) + dum
        END DO
    END DO
    END ASSOCIATE

END SUBROUTINE

!> @brief Subtract damping function term from ESP block
!> @details Compute one-electron Coulomb integrals with the damping function:
!>  \f$ |r-r_C|^{-1} (1 - \beta e^{-\alpha(r-r_C)^2}) \f$
!>  Only the part \f$ - |r-r_C|^{-1} \beta e^{-\alpha(r-r_C)^2}) \f$ is computed here;
!>  the other part is regular Coulomb potential computed elsewhere
!> @param[in]       cp          shell pair data
!> @param[in]       id          current pair of primitives
!> @param[in]       alpha       dumping exponent
!> @param[in]       beta        dumping function scaling factor
!> @param[in]       c           coordinates of the charged particle
!> @param[in]       znuc        particle charge
!> @param[inout]    vblk        block of 1e Coulomb integrals
!
!> @author   Vladimir Mironov
!
!     REVISION HISTORY:
!> @date _Sep, 2018_ Initial release
!
 SUBROUTINE comp_coulomb_dampch_int1_prim(cp,id,alpha,beta,c,znuc,vblk)
!dir$ attributes forceinline :: comp_coulomb_dampch_int1_prim

    TYPE(shpair_t), INTENT(IN) :: cp
    INTEGER, INTENT(IN) :: id
    REAL(REAL64), INTENT(IN) :: alpha, beta, c(3), znuc
    REAL(REAL64), CONTIGUOUS, INTENT(INOUT) :: vblk(:)

    REAL(REAL64) :: xx
    type(rys_root_t) :: ryscomp

    REAL(REAL64) :: dumgij, pcsq, prei, dum, dum1
    INTEGER :: i, j, ij, nx, ny, nz, mx, my, mz, jmax
    real(real64) :: xyzin(0:2*max_ang+1, 0:max_ang+1,3,max_nroots)
!dir$ assume_aligned xyzin : 64
!dir$ assume_aligned vblk : 64

    ASSOCIATE (pp => cp%p(id), &
               iang => cp%iang, jang => cp%jang, &
               inao => cp%inao, jnao => cp%jnao)

    pcsq = pp%aa*sum((pp%r-c)**2)
    xx = pp%aa*pcsq/(pp%aa+alpha)
!   scale DIJ with 1/(aa+alpha) factor
    dumgij = TWOPI/(pp%aa+alpha)
    prei = exp(-(pcsq-xx))

    dum1 = dumgij*pp%expfac*prei*beta

    ryscomp%nroots = cp%nroots
    ryscomp%x = xx
    CALL QGaussRys_damp(ryscomp, cp, id, c, znuc, alpha, xyzin)

    ij = 0
    jmax = jnao
    DO i = 1, inao
        nx = CART_X(i,iang)
        ny = CART_Y(i,iang)
        nz = CART_Z(i,iang)
        IF (cp%iandj) jmax = i
        DO j = 1, jmax
            mx = CART_X(j,jang)
            my = CART_Y(j,jang)
            mz = CART_Z(j,jang)
            ij = ij+1
            dum = sum( xyzin(mx,nx,1,1:cp%nroots)&
                      *xyzin(my,ny,2,1:cp%nroots)&
                      *xyzin(mz,nz,3,1:cp%nroots) )
            vblk(ij) = vblk(ij) - dum1*dum
        END DO
    END DO
    END ASSOCIATE

 END SUBROUTINE

!> @brief Compute sum of 1e Coulomb integrals over primitive shell pair
!> @param[in]       cp          shell pair data
!> @param[in]       id          current pair of primitives
!> @param[in]       c           coordinates of the charged particle
!> @param[in]       den         normalized density matrix block
!> @param[inout]    vsum        sum of Coulomb integrals over pair of primitives
!
!> @author   Vladimir Mironov
!
!     REVISION HISTORY:
!> @date _Oct, 2018_ Initial release
!
 SUBROUTINE comp_coulpot_prim(cp, id, c, den, vsum)
!dir$ attributes inline :: comp_coulpot_prim

    TYPE(shpair_t), INTENT(IN)    :: cp
    INTEGER, INTENT(IN) :: id
    REAL(REAL64),    INTENT(IN)    :: c(3)
    REAL(REAL64),    INTENT(IN)    :: den(:)
    REAL(REAL64),    INTENT(INOUT) :: vsum

    REAL(REAL64) :: xx

    type(rys_root_t) :: ryscomp
    REAL(REAL64) :: tmp
    INTEGER :: i, j, ij, jmax, nx, ny, nz, mx, my, mz
    real(real64) :: xyzin(0:2*max_ang+1, 0:max_ang+1,3,max_nroots)
!dir$ assume_aligned xyzin : 64

    ASSOCIATE (pp => cp%p(id), &
               iang => cp%iang, jang => cp%jang, &
               inao => cp%inao, jnao => cp%jnao)

    xx = pp%aa*sum((pp%r-c)**2)
    ryscomp%nroots = cp%nroots
    ryscomp%x = xx
    CALL QGaussRys(ryscomp, cp, id, c, -1.0d0, xyzin)

    tmp = 0.0
    ij = 0
    jmax = jnao
    DO i = 1, inao
        nx = CART_X(i,iang)
        ny = CART_Y(i,iang)
        nz = CART_Z(i,iang)
        IF (cp%iandj) jmax = i
        DO j = 1, jmax
            mx = CART_X(j,jang)
            my = CART_Y(j,jang)
            mz = CART_Z(j,jang)
            ij = ij+1

            tmp = tmp + den(ij) * sum( xyzin(mx,nx,1,1:cp%nroots)&
                                      *xyzin(my,ny,2,1:cp%nroots)&
                                      *xyzin(mz,nz,3,1:cp%nroots) )

        END DO
    END DO

    vsum = vsum + tmp*pp%expfac*TWOPI*pp%aa1

    END ASSOCIATE

END SUBROUTINE

!> @brief Compute primitive block of 1e Coulomb ESP integrals in FMO method
!> @param[in]       cp          shell pair data
!> @param[in]       id          current pair of primitives
!> @param[inout]    zblk        block of 1e Lz-integrals
!
!> @author   Vladimir Mironov
!
!     REVISION HISTORY:
!> @date _Sep, 2018_ Initial release
!
 SUBROUTINE comp_lz_int1_prim(cp, id, zblk)
!dir$ attributes inline :: int1_lz_prim

    TYPE(shpair_t), INTENT(IN)    :: cp
    INTEGER, INTENT(IN) :: id
    REAL(REAL64), CONTIGUOUS,    INTENT(INOUT) :: zblk(:)

    INTEGER :: i, j, ij, nx, ny, nz, mx, my, mz, jmax

    REAL(REAL64) :: dum2
    real(real64) :: xyzovl(0:max_ang+2,0:max_ang,3)
    real(real64) :: xyzlz(0:max_ang_pad,0:max_ang,2)
!dir$ assume_aligned zblk : 64
!dir$ assume_aligned xyzovl : 64
!dir$ assume_aligned xyzlz : 64

    ASSOCIATE (pp => cp%p(id), &
               iang => cp%iang, jang => cp%jang, &
               inao => cp%inao, jnao => cp%jnao)

    CALL overlap_xyz(cp%ri, cp%rj, pp%r, pp%aa1, iang, jang+1, xyzovl)

    ! j-1
    xyzlz(0,0:iang,1:2) = 0.0
    DO j = 1, jang
        xyzlz(j,0:iang,1:2) = j * xyzovl(j-1,0:iang,1:2)
    END DO

    ij = 0
    jmax = jnao
    DO i = 1, inao
        nx = CART_X(i,iang)
        ny = CART_Y(i,iang)
        nz = CART_Z(i,iang)
        IF (cp%iandj) jmax = i
        DO j = 1, jmax
            mx = CART_X(j,jang)
            my = CART_Y(j,jang)
            mz = CART_Z(j,jang)

            ij = ij+1

            dum2 = xyzovl(mx+1,nx,1)*xyzlz(my,ny,2) - xyzlz(mx,nx,1)*xyzovl(my+1,ny,2)
            zblk(ij) = zblk(ij) + pp%expfac*dum2*xyzovl(mz,nz,3)
        END DO
    END DO
    END ASSOCIATE

END SUBROUTINE

!> @brief Compute primitive block of multipole integrals of order `MOM`
!> @param[in]       cp          shell pair data
!> @param[in]       id          current pair of primitives
!> @param[in]       r           point in space to compute integrals
!> @param[in]       mom         multiplole moment order (1-dipole, 2-quadrupole, 3-octopole)
!> @param[inout]    blk         block of 1e multipole moment integrals
!
!> @author   Vladimir Mironov
!
!     REVISION HISTORY:
!> @date _Sep, 2018_ Initial release
!
 SUBROUTINE comp_mult_int1_prim(cp, id, r, mom, blk)
!dir$ attributes inline :: comp_kin_ovl_int1_prim

    type(shpair_t), intent(in)  :: cp
    integer, intent(in) :: id
    real(real64), contiguous,   intent(in) :: r(:)
    integer, intent(in) :: mom
    real(real64), contiguous,   intent(inout) :: blk(:,:)

    INTEGER :: i, j, nx, ny, nz, mx, my, mz, ij, jmax
    REAL(REAL64) :: tmp(10)
    real(real64) :: xyzmom(3,0:3,0:max_ang,0:max_ang)
!dir$ assume_aligned sblk : 64
!dir$ assume_aligned tblk : 64
!dir$ assume_aligned xyzkin : 64
!dir$ assume_aligned xyzovl : 64

    ASSOCIATE (pp => cp%p(id))
    CALL multipole_xyz(cp%ri, cp%rj, pp%r, pp%aa1, cp%iang, cp%jang, r, mom, xyzmom)

    ij = 0
    jmax = cp%jnao
    DO i = 1, cp%inao
      nx = CART_X(i,cp%iang)
      ny = CART_Y(i,cp%iang)
      nz = CART_Z(i,cp%iang)
      IF (cp%iandj) jmax = i
      DO j = 1, jmax
        mx = CART_X(j,cp%jang)
        my = CART_Y(j,cp%jang)
        mz = CART_Z(j,cp%jang)
        ij = ij+1

        select case (mom)
        case(1)
            tmp(X__) = xyzmom(X__,1,mx,nx)*xyzmom(Y__,0,my,ny)*xyzmom(Z__,0,mz,nz)
            tmp(Y__) = xyzmom(X__,0,mx,nx)*xyzmom(Y__,1,my,ny)*xyzmom(Z__,0,mz,nz)
            tmp(Z__) = xyzmom(X__,0,mx,nx)*xyzmom(Y__,0,my,ny)*xyzmom(Z__,1,mz,nz)
            blk(ij,X__:Z__) =  blk(ij,X__:Z__) + pp%expfac*tmp(X__:Z__)
        case(2)
            tmp(XX_) = xyzmom(X__,2,mx,nx)*xyzmom(Y__,0,my,ny)*xyzmom(Z__,0,mz,nz)
            tmp(YY_) = xyzmom(X__,0,mx,nx)*xyzmom(Y__,2,my,ny)*xyzmom(Z__,0,mz,nz)
            tmp(ZZ_) = xyzmom(X__,0,mx,nx)*xyzmom(Y__,0,my,ny)*xyzmom(Z__,2,mz,nz)
            tmp(XY_) = xyzmom(X__,1,mx,nx)*xyzmom(Y__,1,my,ny)*xyzmom(Z__,0,mz,nz)
            tmp(XZ_) = xyzmom(X__,1,mx,nx)*xyzmom(Y__,0,my,ny)*xyzmom(Z__,1,mz,nz)
            tmp(YZ_) = xyzmom(X__,0,mx,nx)*xyzmom(Y__,1,my,ny)*xyzmom(Z__,1,mz,nz)
            blk(ij,XX_:YZ_) =  blk(ij,XX_:YZ_) + pp%expfac*tmp(XX_:YZ_)
        case(3)
            tmp(XXX) = xyzmom(X__,3,mx,nx)*xyzmom(Y__,0,my,ny)*xyzmom(Z__,0,mz,nz)
            tmp(YYY) = xyzmom(X__,0,mx,nx)*xyzmom(Y__,3,my,ny)*xyzmom(Z__,0,mz,nz)
            tmp(ZZZ) = xyzmom(X__,0,mx,nx)*xyzmom(Y__,0,my,ny)*xyzmom(Z__,3,mz,nz)
            tmp(XXY) = xyzmom(X__,2,mx,nx)*xyzmom(Y__,1,my,ny)*xyzmom(Z__,0,mz,nz)
            tmp(XXZ) = xyzmom(X__,2,mx,nx)*xyzmom(Y__,0,my,ny)*xyzmom(Z__,1,mz,nz)
            tmp(YYX) = xyzmom(X__,1,mx,nx)*xyzmom(Y__,2,my,ny)*xyzmom(Z__,0,mz,nz)
            tmp(YYZ) = xyzmom(X__,0,mx,nx)*xyzmom(Y__,2,my,ny)*xyzmom(Z__,1,mz,nz)
            tmp(ZZX) = xyzmom(X__,1,mx,nx)*xyzmom(Y__,0,my,ny)*xyzmom(Z__,2,mz,nz)
            tmp(ZZY) = xyzmom(X__,0,mx,nx)*xyzmom(Y__,1,my,ny)*xyzmom(Z__,2,mz,nz)
            tmp(XYZ) = xyzmom(X__,1,mx,nx)*xyzmom(Y__,1,my,ny)*xyzmom(Z__,1,mz,nz)
            blk(ij,XXX:XYZ) =  blk(ij,XXX:XYZ) + pp%expfac*tmp(XXX:XYZ)
        case default
            error stop "Max. order of multipole integrals is "// MAX_EL_MOM_S
        end select
      END DO

    END DO

    END ASSOCIATE

 END SUBROUTINE

!> @brief Compute primitive block of multipole
!>        integrals up to an order `MXMOM`
!> @param[in]       cp          shell pair data
!> @param[in]       id          current pair of primitives
!> @param[in]       r           point in space to compute integrals
!> @param[in]       mxmom       multiplole moment order (1-dipole, 2-quadrupole, 3-octopole)
!> @param[inout]    blk         block of 1e multipole moment integrals
!
!> @author   Vladimir Mironov
!
!     REVISION HISTORY:
!> @date _Sep, 2018_ Initial release
!
 SUBROUTINE comp_allmult_int1_prim(cp, id, r, mxmom, blk)
!dir$ attributes inline :: comp_kin_ovl_int1_prim

    type(shpair_t), intent(in)  :: cp
    integer, intent(in) :: id
    real(real64), contiguous,   intent(in) :: r(:)
    integer, intent(in) :: mxmom
    real(real64), contiguous,   intent(inout) :: blk(:,:)

    integer, parameter :: &
        X__ = 1, Y__ = 2, Z__ = 3

    integer, parameter :: &
        XX_ = 4, YY_ = 5, ZZ_ = 6, &
        XY_ = 7, YZ_ = 8, XZ_ = 9

    integer, parameter :: &
        XXX = 10, YYY = 11, ZZZ = 12, &
        XXY = 13, XXZ = 14, &
        YYX = 15, YYZ = 16, &
        ZZX = 17, ZZY = 18, &
        XYZ = 19

    integer, parameter :: BLKDIMS(3) = [3, 3+6, 3+6+10]

    INTEGER :: i, j, nx, ny, nz, mx, my, mz, ij, jmax, blkdim
    REAL(REAL64) :: tmp(19)
    real(real64) :: xyzmom(3,0:3,0:max_ang,0:max_ang)
!dir$ assume_aligned sblk : 64
!dir$ assume_aligned tblk : 64
!dir$ assume_aligned xyzkin : 64
!dir$ assume_aligned xyzovl : 64

    blkdim = blkdims(mxmom)

    ASSOCIATE (pp => cp%p(id))
    CALL multipole_xyz(cp%ri, cp%rj, pp%r, pp%aa1, cp%iang, cp%jang, r, mxmom, xyzmom)

    ij = 0
    jmax = cp%jnao
    DO i = 1, cp%inao
      nx = CART_X(i,cp%iang)
      ny = CART_Y(i,cp%iang)
      nz = CART_Z(i,cp%iang)
      IF (cp%iandj) jmax = i
      DO j = 1, jmax
        mx = CART_X(j,cp%jang)
        my = CART_Y(j,cp%jang)
        mz = CART_Z(j,cp%jang)
        ij = ij+1

        tmp(X__) = xyzmom(X__,1,mx,nx)*xyzmom(Y__,0,my,ny)*xyzmom(Z__,0,mz,nz)
        tmp(Y__) = xyzmom(X__,0,mx,nx)*xyzmom(Y__,1,my,ny)*xyzmom(Z__,0,mz,nz)
        tmp(Z__) = xyzmom(X__,0,mx,nx)*xyzmom(Y__,0,my,ny)*xyzmom(Z__,1,mz,nz)
        if (mxmom>1) then
            tmp(XX_) = xyzmom(X__,2,mx,nx)*xyzmom(Y__,0,my,ny)*xyzmom(Z__,0,mz,nz)
            tmp(YY_) = xyzmom(X__,0,mx,nx)*xyzmom(Y__,2,my,ny)*xyzmom(Z__,0,mz,nz)
            tmp(ZZ_) = xyzmom(X__,0,mx,nx)*xyzmom(Y__,0,my,ny)*xyzmom(Z__,2,mz,nz)
            tmp(XY_) = xyzmom(X__,1,mx,nx)*xyzmom(Y__,1,my,ny)*xyzmom(Z__,0,mz,nz)
            tmp(YZ_) = xyzmom(X__,0,mx,nx)*xyzmom(Y__,1,my,ny)*xyzmom(Z__,1,mz,nz)
            tmp(XZ_) = xyzmom(X__,1,mx,nx)*xyzmom(Y__,0,my,ny)*xyzmom(Z__,1,mz,nz)
        end if
        if (mxmom>2) then
            tmp(XXX) = xyzmom(X__,3,mx,nx)*xyzmom(Y__,0,my,ny)*xyzmom(Z__,0,mz,nz)
            tmp(YYY) = xyzmom(X__,0,mx,nx)*xyzmom(Y__,3,my,ny)*xyzmom(Z__,0,mz,nz)
            tmp(ZZZ) = xyzmom(X__,0,mx,nx)*xyzmom(Y__,0,my,ny)*xyzmom(Z__,3,mz,nz)
            tmp(XXY) = xyzmom(X__,2,mx,nx)*xyzmom(Y__,1,my,ny)*xyzmom(Z__,0,mz,nz)
            tmp(XXZ) = xyzmom(X__,2,mx,nx)*xyzmom(Y__,0,my,ny)*xyzmom(Z__,1,mz,nz)
            tmp(YYX) = xyzmom(X__,1,mx,nx)*xyzmom(Y__,2,my,ny)*xyzmom(Z__,0,mz,nz)
            tmp(YYZ) = xyzmom(X__,0,mx,nx)*xyzmom(Y__,2,my,ny)*xyzmom(Z__,1,mz,nz)
            tmp(ZZX) = xyzmom(X__,1,mx,nx)*xyzmom(Y__,0,my,ny)*xyzmom(Z__,2,mz,nz)
            tmp(ZZY) = xyzmom(X__,0,mx,nx)*xyzmom(Y__,1,my,ny)*xyzmom(Z__,2,mz,nz)
            tmp(XYZ) = xyzmom(X__,1,mx,nx)*xyzmom(Y__,1,my,ny)*xyzmom(Z__,1,mz,nz)
        end if
        blk(ij,1:blkdim) =  blk(ij,1:blkdim) + pp%expfac*tmp(1:blkdim)
      END DO

    END DO

    END ASSOCIATE

 END SUBROUTINE

!--------------------------------------------------------------------------------
!       ONE-ELECTRON DERIVATIVES CALCULATION (PRIMITIVE GAUSSIANS)
!--------------------------------------------------------------------------------

!> @brief Compute 1e overlap contribution to the gradient
!> @param[in]       cp          shell pair data
!> @param[in]       dij         density matrix block
!> @param[inout]    de          dimension(3), contribution to gradient
!
!> @author   Vladimir Mironov
!
!     REVISION HISTORY:
!> @date _Sep, 2018_ Initial release
!
 SUBROUTINE comp_overlap_der1(cp, dij, de)
!dir$ attributes inline :: comp_overlap_der1
    TYPE(shpair_t), INTENT(IN) :: cp
    REAL(REAL64), INTENT(IN) :: dij(:,:)
    REAL(REAL64), CONTIGUOUS, INTENT(INOUT) :: de(:)

    REAL(REAL64) :: der(3), de_loc(3)
    INTEGER :: i, j, k, ix, iy, iz, jx, jy, jz
    real(real64) :: ovl_int(0:max_ang,0:max_ang+3,3)
    real(real64) :: ovl_der(0:max_ang,0:max_ang,3)

    DO k = 1, cp%numpairs
    ASSOCIATE (pp => cp%p(k), &
               iang => cp%iang, jang => cp%jang, &
               inao => cp%inao, jnao => cp%jnao)
    ! compute overlap [i+1|j]
    CALL overlap_xyz(cp%ri, cp%rj, pp%r, pp%aa1, iang+1, jang, ovl_int)

    ! compute 1D overlap derivatives [i|j]
    CALL der_kinovl_xyz(ovl_der,ovl_int,iang,jang,pp%ai)

    ! assemble overlap contribution to the gradient
    de_loc = 0.0
    DO i = 1, inao
        ix = CART_X(i,iang)
        iy = CART_Y(i,iang)
        iz = CART_Z(i,iang)
        DO j = 1, jnao
            jx = CART_X(j,jang)
            jy = CART_Y(j,jang)
            jz = CART_Z(j,jang)
            der(1) = ovl_der(jx,ix,1) * ovl_int(jy,iy,2) * ovl_int(jz,iz,3)
            der(2) = ovl_int(jx,ix,1) * ovl_der(jy,iy,2) * ovl_int(jz,iz,3)
            der(3) = ovl_int(jx,ix,1) * ovl_int(jy,iy,2) * ovl_der(jz,iz,3)
            de_loc = de_loc + der*dij(i,j)
        END DO
    END DO
    ! add scaled contribution to gradient
    de = de + de_loc*pp%expfac
    END ASSOCIATE
    END DO
 END SUBROUTINE

!> @brief Compute 1e kinetic contribution to the gradient
!> @param[in]       cp          shell pair data
!> @param[in]       dij         density matrix block
!> @param[inout]    de          dimension(3), contribution to gradient
!
!> @author   Vladimir Mironov
!
!     REVISION HISTORY:
!> @date _Sep, 2018_ Initial release
!
 SUBROUTINE comp_kinetic_der1(cp, dij, de)
!dir$ attributes inline :: comp_kinetic_der1
    TYPE(shpair_t), INTENT(IN) :: cp
    REAL(REAL64), INTENT(IN) :: dij(:,:)
    REAL(REAL64), CONTIGUOUS, INTENT(INOUT) :: de(:)

    REAL(REAL64) :: der(3), de_loc(3)
    INTEGER :: i, j, k, ix, iy, iz, jx, jy, jz
    real(real64) :: ovl_int(0:max_ang,0:max_ang+3,3)
    real(real64) :: ovl_der(0:max_ang,0:max_ang,3)

    real(real64) :: kin_int(0:max_ang,0:max_ang+1,3)
    real(real64) :: kin_der(0:max_ang,0:max_ang,3)


    DO k = 1, cp%numpairs
    ASSOCIATE (pp => cp%p(k), &
               iang => cp%iang, jang => cp%jang, &
               inao => cp%inao, jnao => cp%jnao)
    ! compute overlap [i+3|j]
    CALL overlap_xyz(cp%ri, cp%rj, pp%r, pp%aa1, iang+3, jang, ovl_int)

    ! compute 1D kinetic [i+1|j]
    CALL kinetic_xyz_i(kin_int,ovl_int,iang+1,jang,pp%ai)

    ! compute 1D overlap derivatives [i|j]
    CALL der_kinovl_xyz(ovl_der,ovl_int,iang,jang,pp%ai)

    ! compute 1D kinetic derivatives [i|j]
    CALL der_kinovl_xyz(kin_der,kin_int,iang,jang,pp%ai)

    ! assemble 3D K.E. derivatives from 1D integrals and derivatives
    de_loc = 0.0
    DO i = 1, inao
        ix = CART_X(i,iang)
        iy = CART_Y(i,iang)
        iz = CART_Z(i,iang)
        DO j = 1, jnao
            jx = CART_X(j,jang)
            jy = CART_Y(j,jang)
            jz = CART_Z(j,jang)
            der(1) = kin_der(jx,ix,1) * ovl_int(jy,iy,2) * ovl_int(jz,iz,3) + &
                     ovl_der(jx,ix,1) * kin_int(jy,iy,2) * ovl_int(jz,iz,3) + &
                     ovl_der(jx,ix,1) * ovl_int(jy,iy,2) * kin_int(jz,iz,3)
            der(2) = kin_int(jx,ix,1) * ovl_der(jy,iy,2) * ovl_int(jz,iz,3) + &
                     ovl_int(jx,ix,1) * kin_der(jy,iy,2) * ovl_int(jz,iz,3) + &
                     ovl_int(jx,ix,1) * ovl_der(jy,iy,2) * kin_int(jz,iz,3)
            der(3) = kin_int(jx,ix,1) * ovl_int(jy,iy,2) * ovl_der(jz,iz,3) + &
                     ovl_int(jx,ix,1) * kin_int(jy,iy,2) * ovl_der(jz,iz,3) + &
                     ovl_int(jx,ix,1) * ovl_int(jy,iy,2) * kin_der(jz,iz,3)
            de_loc = de_loc + der*dij(i,j)
        END DO
    END DO
    ! add scaled contribution to gradient
    de = de + de_loc*pp%expfac
    END ASSOCIATE
    END DO
 END SUBROUTINE

!> @brief Compute 1e Coulomb contribution to the gradient (v.r.t. shifts of
!>  shell's centers)
!> @param[in]       nroots      roots for GaussRys
!> @param[in]       cp          shell pair data
!> @param[in]       c           coordinates of the charged particle
!> @param[in]       znuc        particle charge
!> @param[in]       dij         density matrix block
!> @param[inout]    dernuc      dimension(3), contribution to gradient
!
!> @author   Vladimir Mironov
!
!     REVISION HISTORY:
!> @date _Sep, 2018_ Initial release
!
 SUBROUTINE comp_coulomb_der1(cp, c, znuc, dij, dernuc)
!dir$ attributes inline :: comp_coulomb_der1
    TYPE(shpair_t), INTENT(IN) :: cp
    REAL(REAL64), INTENT(IN) :: c(3), znuc
    REAL(REAL64), INTENT(IN) :: dij(:,:)
    REAL(REAL64), INTENT(OUT) :: dernuc(3)

    REAL(REAL64) :: xx
    type(rys_root_t) :: ryscomp
    REAL(REAL64) :: der(3), fac, detmp(3)
    INTEGER :: id, i, j, ix, iy, iz, jx, jy, jz
    real(real64) :: xyzin(0:2*max_ang+1, 0:max_ang+1,3,max_nroots)
    real(real64) :: dxyzc(0:max_ang_pad,0:max_ang,3,max_nroots)
!dir$ assume_aligned xyzin : 64
!dir$ assume_aligned dxyzc : 64

    dernuc = 0.0

    DO id = 1, cp%numpairs

        ASSOCIATE (pp => cp%p(id), &
                   iang => cp%iang, jang => cp%jang, &
                   inao => cp%inao, jnao => cp%jnao)

        xx = pp%aa*sum((pp%r-c)**2)

        ryscomp%nroots = cp%nroots
        ryscomp%x = xx
        CALL QGaussRys(ryscomp, cp, id, c, znuc, xyzin, 1)

        CALL der_coul_xyz(dxyzc,xyzin,iang,jang,pp%ai,cp%nroots)

        fac = pp%expfac*TWOPI*pp%aa1

        detmp = 0.0

        DO i = 1, inao
            ix = CART_X(i,iang)
            iy = CART_Y(i,iang)
            iz = CART_Z(i,iang)
            DO j = 1, jnao
                jx = CART_X(j,jang)
                jy = CART_Y(j,jang)
                jz = CART_Z(j,jang)

                der(1) = sum( dxyzc(jx,ix,1,1:cp%nroots)&
                             *xyzin(jy,iy,2,1:cp%nroots)&
                             *xyzin(jz,iz,3,1:cp%nroots) )

                der(2) = sum( xyzin(jx,ix,1,1:cp%nroots)&
                             *dxyzc(jy,iy,2,1:cp%nroots)&
                             *xyzin(jz,iz,3,1:cp%nroots) )

                der(3) = sum( xyzin(jx,ix,1,1:cp%nroots)&
                             *xyzin(jy,iy,2,1:cp%nroots)&
                             *dxyzc(jz,iz,3,1:cp%nroots) )

                detmp = detmp + der*dij(i,j)

             END DO
        END DO
        dernuc = dernuc + detmp*fac
        END ASSOCIATE
    END DO

 END SUBROUTINE

!> @brief Compute 1e Hellmann-Feynman contribution to the gradient
!> @param[in]       nroots      roots for GaussRys
!> @param[in]       cp          shell pair data
!> @param[in]       c           coordinates of the charged particle
!> @param[in]       znuc        particle charge
!> @param[in]       dij         density matrix block
!> @param[inout]    derhf       dimension(3), contribution to gradient
!
!> @author   Vladimir Mironov
!
!     REVISION HISTORY:
!> @date _Sep, 2018_ Initial release
!
 SUBROUTINE comp_coulomb_helfeyder1(cp, c, znuc, dij, derhf)
!dir$ attributes inline :: comp_coulomb_helfeyder1
    TYPE(shpair_t), INTENT(IN) :: cp
    REAL(REAL64), INTENT(IN) :: c(3), znuc
    REAL(REAL64), INTENT(IN) :: dij(:,:)
    REAL(REAL64), CONTIGUOUS, INTENT(OUT) :: derhf(:)

    REAL(REAL64) :: xx
    type(rys_root_t) :: ryscomp
    REAL(REAL64) :: ric(3)
    REAL(REAL64) :: der(3), fac
    INTEGER :: id, i, j, ix, iy, iz, jx, jy, jz
    real(real64) :: xyzin(0:2*max_ang+1, 0:max_ang+1,3,max_nroots)
    real(real64) :: dxyzc(0:max_ang_pad,0:max_ang,3,max_nroots)
!dir$ assume_aligned xyzin : 64
!dir$ assume_aligned dxyzc : 64

    derhf = 0.0

    ric = cp%ri(:3) - c(:3)

    DO id = 1, cp%numpairs

        ASSOCIATE (pp => cp%p(id), &
                   iang => cp%iang, jang => cp%jang, &
                   inao => cp%inao, jnao => cp%jnao)

        xx = pp%aa*sum((pp%r-c)**2)

        fac = pp%expfac*TWOPI*2

        ryscomp%nroots = cp%nroots
        ryscomp%x = xx
        CALL DQGaussRys(ryscomp, cp, id, c, znuc, xyzin)

        CALL der_helfey_xyz(dxyzc,xyzin,iang,jang,ric,cp%nroots)

        DO i = 1, inao
            ix = CART_X(i,iang)
            iy = CART_Y(i,iang)
            iz = CART_Z(i,iang)
            DO j = 1, jnao
                jx = CART_X(j,jang)
                jy = CART_Y(j,jang)
                jz = CART_Z(j,jang)

                der(1) = sum(dxyzc(jx,ix,1,1:cp%nroots)&
                            *xyzin(jy,iy,2,1:cp%nroots)&
                            *xyzin(jz,iz,3,1:cp%nroots))

                der(2) = sum(xyzin(jx,ix,1,1:cp%nroots)&
                            *dxyzc(jy,iy,2,1:cp%nroots)&
                            *xyzin(jz,iz,3,1:cp%nroots))

                der(3) = sum(xyzin(jx,ix,1,1:cp%nroots)&
                            *xyzin(jy,iy,2,1:cp%nroots)&
                            *dxyzc(jz,iz,3,1:cp%nroots))

                derhf = derhf + der*dij(i,j)*fac
             END DO
        END DO
        END ASSOCIATE
    END DO

 END SUBROUTINE

!> @brief Compute 1e Ewald long-range contribution to the gradient (v.r.t. shifts of
!>  shell's centers)
!> @param[in]       nroots      roots for GaussRys
!> @param[in]       cp          shell pair data
!> @param[in]       c           coordinates of the charged particle
!> @param[in]       znuc        particle charge
!> @param[in]       dij         density matrix block
!> @param[in]       omega       Ewald splitting parameter
!> @param[inout]    dernuc      dimension(3), contribution to gradient
!
!> @author   Vladimir Mironov
!
!     REVISION HISTORY:
!> @date _Sep, 2018_ Initial release
!
 SUBROUTINE comp_ewaldlr_der1(cp, c, znuc, dij, omega, dernuc)
!dir$ attributes inline :: comp_ewaldlr_der1
    TYPE(shpair_t), INTENT(IN) :: cp
    REAL(REAL64), INTENT(IN) :: c(3), znuc
    REAL(REAL64), INTENT(IN) :: dij(:,:)
    REAL(REAL64), INTENT(IN)   :: omega
    REAL(REAL64), INTENT(OUT) :: dernuc(3)

    REAL(REAL64) :: xx
    type(rys_root_t) :: ryscomp

    REAL(REAL64) :: der(3), fac, detmp(3), xfac
    INTEGER :: id, i, j, ix, iy, iz, jx, jy, jz
    real(real64) :: xyzin(0:2*max_ang+1, 0:max_ang+1,3,max_nroots)
    real(real64) :: dxyzc(0:max_ang_pad,0:max_ang,3,max_nroots)
!dir$ assume_aligned xyzin : 64
!dir$ assume_aligned dxyzc : 64

    !dernuc = 0.0

    DO id = 1, cp%numpairs

        ASSOCIATE (pp => cp%p(id), &
                   iang => cp%iang, jang => cp%jang, &
                   inao => cp%inao, jnao => cp%jnao)

        xfac = omega*omega/(pp%aa+omega*omega)
        xx = pp%aa* sum((pp%r-c)**2) * xfac
        ryscomp%nroots = cp%nroots
        ryscomp%x = xx
        CALL QGaussRysEw(ryscomp, cp, id, c, znuc, xfac, xyzin, 1)

        CALL der_coul_xyz(dxyzc,xyzin,iang,jang,pp%ai,cp%nroots)

        fac = pp%expfac*TWOPI*pp%aa1 * sqrt(xfac)

        detmp = 0.0

        DO i = 1, inao
            ix = CART_X(i,iang)
            iy = CART_Y(i,iang)
            iz = CART_Z(i,iang)
            DO j = 1, jnao
                jx = CART_X(j,jang)
                jy = CART_Y(j,jang)
                jz = CART_Z(j,jang)

                der(1) = sum( dxyzc(jx,ix,1,1:cp%nroots)&
                             *xyzin(jy,iy,2,1:cp%nroots)&
                             *xyzin(jz,iz,3,1:cp%nroots) )

                der(2) = sum( xyzin(jx,ix,1,1:cp%nroots)&
                             *dxyzc(jy,iy,2,1:cp%nroots)&
                             *xyzin(jz,iz,3,1:cp%nroots) )

                der(3) = sum( xyzin(jx,ix,1,1:cp%nroots)&
                             *xyzin(jy,iy,2,1:cp%nroots)&
                             *dxyzc(jz,iz,3,1:cp%nroots) )

                detmp = detmp + der*dij(i,j)

             END DO
        END DO
        dernuc = dernuc + detmp*fac
        END ASSOCIATE
    END DO

 END SUBROUTINE

!> @brief Compute Ewald long-range 1e Hellmann-Feynman contribution
!>  to the gradient
!> @param[in]       nroots      roots for GaussRys
!> @param[in]       cp          shell pair data
!> @param[in]       c           coordinates of the charged particle
!> @param[in]       znuc        particle charge
!> @param[in]       dij         density matrix block
!> @param[in]       omega       Ewald splitting parameter
!> @param[inout]    derhf       dimension(3), contribution to gradient
!
!> @author   Vladimir Mironov
!
!     REVISION HISTORY:
!> @date _Sep, 2018_ Initial release
!
 SUBROUTINE comp_ewaldlr_helfeyder1(cp, c, znuc, dij, omega, derhf)
!dir$ attributes inline :: comp_ewaldlr_helfeyder1
    TYPE(shpair_t), INTENT(IN) :: cp
    REAL(REAL64), INTENT(IN) :: c(3), znuc
    REAL(REAL64), INTENT(IN) :: dij(:,:)
    REAL(REAL64), INTENT(IN)   :: omega
    REAL(REAL64), CONTIGUOUS, INTENT(OUT) :: derhf(:)

    REAL(REAL64) :: xx
    type(rys_root_t) :: ryscomp

    REAL(REAL64) :: ric(3)
    REAL(REAL64) :: der(3), fac, xfac
    INTEGER :: id, i, j, ix, iy, iz, jx, jy, jz
    real(real64) :: xyzin(0:2*max_ang+1, 0:max_ang+1,3,max_nroots)
    real(real64) :: dxyzc(0:max_ang_pad,0:max_ang,3,max_nroots)
!dir$ assume_aligned xyzin : 64
!dir$ assume_aligned dxyzc : 64

    derhf = 0.0

    ric = cp%ri(:3) - c(:3)

    DO id = 1, cp%numpairs

        ASSOCIATE (pp => cp%p(id), &
                   iang => cp%iang, jang => cp%jang, &
                   inao => cp%inao, jnao => cp%jnao)

        xfac = omega*omega/(pp%aa+omega*omega)
        xx = pp%aa* sum((pp%r-c)**2) * xfac
        ryscomp%nroots = cp%nroots
        ryscomp%x = xx
        CALL DQGaussRysEw(ryscomp, cp, id, c, znuc, xfac, xyzin)

        CALL der_helfey_xyz(dxyzc,xyzin,iang,jang,ric,cp%nroots)

        fac = pp%expfac*TWOPI*2 * sqrt(xfac)

        DO i = 1, inao
            ix = CART_X(i,iang)
            iy = CART_Y(i,iang)
            iz = CART_Z(i,iang)
            DO j = 1, jnao
                jx = CART_X(j,jang)
                jy = CART_Y(j,jang)
                jz = CART_Z(j,jang)

                der(1) = sum(dxyzc(jx,ix,1,1:cp%nroots)&
                            *xyzin(jy,iy,2,1:cp%nroots)&
                            *xyzin(jz,iz,3,1:cp%nroots))

                der(2) = sum(xyzin(jx,ix,1,1:cp%nroots)&
                            *dxyzc(jy,iy,2,1:cp%nroots)&
                            *xyzin(jz,iz,3,1:cp%nroots))

                der(3) = sum(xyzin(jx,ix,1,1:cp%nroots)&
                            *xyzin(jy,iy,2,1:cp%nroots)&
                            *dxyzc(jz,iz,3,1:cp%nroots))

                derhf = derhf + der*dij(i,j)*fac
             END DO
        END DO
        END ASSOCIATE
    END DO

 END SUBROUTINE

!--------------------------------------------------------------------------------
!       K.E. AND OVERLAP 1D INTEGRALS
!--------------------------------------------------------------------------------

!> @brief Compute 1D overlap integrals
!> @details Return block of 1D integrals, dimensions: (Lj,Li,XYZ)
!> @param[in]       ri          coordinates of first shell center
!> @param[in]       rj          coordinates of second shell center
!> @param[in]       rij         coordinates of shell-pair center of charge
!> @param[in]       aa1         inverse total exponent
!> @param[in]       li          max angular momentum for the first shell center
!> @param[in]       lj          max angular momentum for the second shell center
!> @param[inout]    xyzovl      block of 1D overlap integrals
!
!> @author   Vladimir Mironov
!
!     REVISION HISTORY:
!> @date _Sep, 2018_ Initial release
!
 SUBROUTINE overlap_xyz(ri, rj, rij, aa1, li, lj, xyzovl)
!dir$ attributes forceinline :: overlap_xyz

    REAL(REAL64), CONTIGUOUS, INTENT(IN)   :: ri(:), rj(:), rij(:)
    REAL(REAL64), INTENT(IN)   :: aa1
    INTEGER, INTENT(IN)   :: li, lj
    REAL(REAL64), CONTIGUOUS, INTENT(OUT)  :: xyzovl(0:,0:,:)

    INTEGER :: i, j
    REAL(REAL64) :: taa, oint(3)

    taa = sqrt(aa1)

    DO i = 0, li
        DO j = 0, lj
            CALL doQuadGaussHermite(oint, taa, rij(:3), &
                                    ri(:3), rj(:3), i, j)
            xyzovl(j,i,:) = oint*taa
        END DO
    END DO

 END SUBROUTINE

!> @brief Kinetic energy integrals, recursion over first shell
!> @details Compute K.E.I. from overlap integrals using recurrence
!>  over first shell
!> @param[out] xyzt 1D kinetic energy integrals
!> @param[in]  xyzs 1D overlap integrals
!> @param[in]  ni   number of points for the 1st shell quad.
!> @param[in]  nj   number of points for the 2nd shell quad.
!> @param[in]  ai   first shell exponent
!> @author   Vladimir Mironov
!
!> @note Before running this routine, first you need to
!>  compute overlap integrals for angular momentums (Li+2, Lj)
!
!     REVISION HISTORY:
!> @date _Sep, 2018_ Initial release
!
 SUBROUTINE kinetic_xyz_i(xyzt,xyzs,ni,nj,ai)
!dir$ attributes forceinline :: kinetic_xyz_i
    REAL(REAL64), CONTIGUOUS, INTENT(OUT) :: xyzt(0:,0:,:)
    REAL(REAL64), CONTIGUOUS, INTENT(IN)  :: xyzs(0:,0:,:)
    REAL(REAL64), INTENT(IN)  :: ai
    INTEGER, INTENT(IN) :: ni, nj

    INTEGER :: i
    REAL(REAL64) :: fact1, fact2

    xyzt(0:nj,0,:) = (xyzs(0:nj,0,:) - 2*ai*xyzs(0:nj,2,:))*ai

    IF (ni==0) RETURN

    xyzt(0:nj,1,:) = (xyzs(0:nj,1,:)*3.0D0 - 2*ai*xyzs(0:nj,3,:))*ai

    IF (ni==1) RETURN

    DO i = 2, ni
      fact1 = 2*i+1
      fact2 = real(i*(i-1)/2,REAL64)
      xyzt(0:nj,i,:) = (xyzs(0:nj,i,:)*fact1 - 2*ai*xyzs(0:nj,i+2,:))*ai - xyzs(0:nj,i-2,:)*fact2
    END DO

 END SUBROUTINE

!> @brief Kinetic energy integrals, recursion over second shell
!> @details Compute K.E.I. from overlap integrals using recurrence
!>  over second shell
!> @param[out] xyzt 1D kinetic energy integrals
!> @param[in]  xyzs 1D overlap integrals
!> @param[in]  ni   number of points for the 1st shell quad.
!> @param[in]  nj   number of points for the 2nd shell quad.
!> @param[in]  aj   second shell exponent
!> @author   Vladimir Mironov
!
!> @note Before running this routine, first you need to
!>  compute overlap integrals for angular momentums (Li+2, Lj)
!
!     REVISION HISTORY:
!> @date _Sep, 2018_ Initial release
!
 SUBROUTINE kinetic_xyz_j(xyzt,xyzs,ni,nj,aj)
!dir$ attributes forceinline :: kinetic_xyz_j
    REAL(REAL64), CONTIGUOUS, INTENT(OUT) :: xyzt(0:,0:,:)
    REAL(REAL64), CONTIGUOUS, INTENT(IN)  :: xyzs(0:,0:,:)
    REAL(REAL64), INTENT(IN)  :: aj
    INTEGER, INTENT(IN) :: ni, nj

    INTEGER :: j
    REAL(REAL64) :: fact1, fact2

    xyzt(0,0:ni,:) = (xyzs(0,0:ni,:) - 2*aj*xyzs(2,0:ni,:))*aj

    IF (nj==0) RETURN

    xyzt(1,0:ni,:) = (xyzs(1,0:ni,:)*3.0D0 - 2*aj*xyzs(3,0:ni,:))*aj

    IF (nj==1) RETURN

    DO j = 2, nj
      fact1 = 2*j+1
      fact2 = real(j*(j-1)/2,REAL64)
      xyzt(j,0:ni,:) = (xyzs(j,0:ni,:)*fact1 - 2*aj*xyzs(j+2,0:ni,:))*aj - xyzs(j-2,0:ni,:)*fact2
    END DO

 END SUBROUTINE

!--------------------------------------------------------------------------------
!       Multipole moment integrals
!--------------------------------------------------------------------------------

!> @brief Compute multipole moment integrals using Gauss-Hermite quadrature
!> @details Return block of 1D integrals, dimensions: (Lj,Li,XYZ)
!> @param[in]       ri          coordinates of first shell center
!> @param[in]       rj          coordinates of second shell center
!> @param[in]       rij         coordinates of shell-pair center of charge
!> @param[in]       aa1         inverse total exponent
!> @param[in]       li          max angular momentum for the first shell center
!> @param[in]       lj          max angular momentum for the second shell center
!> @param[in]       r           origin of multipole moment integrals
!> @param[in]       mxmom       max order of multipole moment integrals
!> @param[inout]    xyzints     block of mutipole moments integrals
!
!> @author   Vladimir Mironov
!
!     REVISION HISTORY:
!> @date _Sep, 2018_ Initial release
!
 SUBROUTINE multipole_xyz(ri, rj, rij, aa1, li, lj, r, mxmom, xyzints)
!dir$ attributes forceinline :: overlap_xyz

    REAL(REAL64), CONTIGUOUS, INTENT(IN)   :: ri(:), rj(:), rij(:), r(:)
    REAL(REAL64), INTENT(IN)   :: aa1
    INTEGER, INTENT(IN)   :: li, lj
    integer, intent(in) :: mxmom
    REAL(REAL64), CONTIGUOUS, INTENT(OUT)  :: xyzints(:,0:,0:,0:)

    INTEGER :: i, j
    REAL(REAL64) :: taa, ppint(3,0:MAX_EL_MOM)

    taa = sqrt(aa1)

    DO i = 0, li
        DO j = 0, lj
            CALL mulQuadGaussHermite(ppint, taa, rij(:3), &
                                    ri(:3), rj(:3), r(:3), i, j, mxmom)
            xyzints(:,0:mxmom,j,i) = ppint(:,0:mxmom)*taa

        END DO
    END DO

 END SUBROUTINE

!--------------------------------------------------------------------------------
!       DERIVATIVE CODE 1D INTEGRALS
!--------------------------------------------------------------------------------

!> @brief Compute derivatives of 1D Coulomb integrals v.r.t. shifts of shell centers
!> @details Derivatives are computed using following equation:
!>  \f$ D(L_i,L_j) = 2\alpha_i I(L_i+1,L_j) - L_i I(L_i-1,L_j) \f$
!> @param[out] dxyzdi    1D Coulomb integral derivatives (dims: (Lj,Li,XYZ,NRoots)
!> @param[in]  xyzin     1D Coulomb integrals (dims: (Lj,Li,XYZ,NRoots)
!> @param[in]  lit       angular momentum of the 1st shell + 1
!> @param[in]  ljt       angular momentum of the 2nd shell + 1
!> @param[in]  ai        exponent of the first shell
!> @param[in]  nroots    number of roots in Gauss-Rys quadrature
!> @author   Vladimir Mironov
!
!> @note based on DERI from grd1.src
!
!     REVISION HISTORY:
!> @date _Sep, 2018_ Initial release
!
 SUBROUTINE der_coul_xyz(dxyzdi,xyzin,lit,ljt,ai,nroots)
!dir$ attributes forceinline :: der_coul_xyz
    REAL(REAL64), INTENT(IN) ::  ai
    REAL(REAL64), CONTIGUOUS, INTENT(IN) ::  xyzin(0:,0:,:,:)
    REAL(REAL64), CONTIGUOUS, INTENT(OUT) :: dxyzdi(0:,0:,:,:)
    INTEGER, INTENT(IN) :: lit, ljt, nroots

    INTEGER :: i
!dir$ assume_aligned xyzin : 64
!dir$ assume_aligned dxyzdi : 64

    dxyzdi(0:ljt,0:lit,1:3,1:nroots) = 2*ai * xyzin(0:ljt,1:lit+1,1:3,1:nroots)

    DO i = 1, lit
        dxyzdi(0:ljt,i,1:3,1:nroots) = dxyzdi(0:ljt,i,1:3,1:nroots) - i*xyzin(0:ljt,i-1,1:3,1:nroots)
    END DO

 END SUBROUTINE

!> @brief Compute derivatives of 1D Coulomb integrals v.r.t. shifts of the nuclei
!>  (Hellman-Feynman term)
!>  Derivatives are computed using following equation:
!>  \f$ D(L_i,L_j) = I(L_i+1,L_j) - (r_i - r_c) I(L_i,L_j) \f$
!> @param[out] dxyzdc    1D Coulomb integral derivatives (dims: (Lj,Li,XYZ,NRoots)
!> @param[in]  xyzin     1D Coulomb integrals (dims: (Lj,Li,XYZ,NRoots)
!> @param[in]  lit       angular momentum of the 1st shell + 1
!> @param[in]  ljt       angular momentum of the 2nd shell + 1
!> @param[in]  ric       (Ri-Rij)xyz
!> @param[in]  nroots    number of roots in Gauss-Rys quadrature
!> @author   Vladimir Mironov
!
!     REVISION HISTORY:
!> @date _Sep, 2018_ Initial release
!
 SUBROUTINE der_helfey_xyz(dxyzdc,xyzin,lit,ljt,ric,nroots)
!dir$ attributes forceinline :: der_helfey_xyz
    REAL(REAL64), INTENT(IN) ::  ric(3)
    REAL(REAL64), CONTIGUOUS, INTENT(IN) ::  xyzin(0:,0:,:,:)
    REAL(REAL64), CONTIGUOUS, INTENT(OUT) :: dxyzdc(0:,0:,:,:)
    INTEGER, INTENT(IN) :: lit, ljt, nroots
!dir$ assume_aligned xyzin : 64

    dxyzdc(0:ljt,0:lit,1,1:nroots) = xyzin(0:ljt,1:lit+1,1,1:nroots) + ric(1)*xyzin(0:ljt,0:lit,1,1:nroots)
    dxyzdc(0:ljt,0:lit,2,1:nroots) = xyzin(0:ljt,1:lit+1,2,1:nroots) + ric(2)*xyzin(0:ljt,0:lit,2,1:nroots)
    dxyzdc(0:ljt,0:lit,3,1:nroots) = xyzin(0:ljt,1:lit+1,3,1:nroots) + ric(3)*xyzin(0:ljt,0:lit,3,1:nroots)

 END SUBROUTINE

!> @brief 1e overlap and kinetic energy integrals differentiation
!> @param[out] dxyz 1D derivatives
!> @param[in]  xyz  1D integrals
!> @param[in]  lit  angular momentum of the 1st shell + 1
!> @param[in]  ljt  angular momentum of the 2nd shell + 1
!> @param[in]  ai   exponent of the first shell
!> @author   Vladimir Mironov
!
!     REVISION HISTORY:
!> @date _Sep, 2018_ Initial release
!
 SUBROUTINE der_kinovl_xyz(dxyz,xyz,lit,ljt,ai)
!dir$ attributes forceinline :: der_kinovl_xyz
    REAL(REAL64), INTENT(IN) ::  ai
    REAL(REAL64), CONTIGUOUS, INTENT(IN) ::  xyz(0:,0:,:)
    REAL(REAL64), CONTIGUOUS, INTENT(OUT) :: dxyz(0:,0:,:)
    INTEGER, INTENT(IN) :: lit, ljt
    INTEGER :: i

    dxyz(0:ljt,0:lit,:) =  2*ai * xyz(0:ljt,1:lit+1,:)

    !IF (lit==1) RETURN

    DO i = 1, lit
        dxyz(0:ljt,i,:) = dxyz(0:ljt,i,:) - i*xyz(0:ljt,i-1,:)
    END DO

 END SUBROUTINE

!--------------------------------------------------------------------------------
!       GAUSS-RYS QUADRATURE RELATED ROUTINES
!--------------------------------------------------------------------------------

!> @brief Compute 1D integrals for 1e Coulomb integrals
!> @details In this implementation 1D integrals at Rys abscissae are
!>  computed using VRR and HRR recurrences
!> @note The common factor for the integral block is \f$ 2\Pi \f$
!> @param[in]  nroots   roots for GaussRys
!> @param[in]  cp       shell pair data
!> @param[in]  id       current pair of primitives
!> @param[in]  c        coordinates of the charged particle
!> @param[in]  znuc     charge of the particle
!> @param[out] xyzin    array of 1D integrals
!> @param[in]  igrd     [opt] flag indicating that integral derivatives are needed
!> @author   Vladimir Mironov
!
!     REVISION HISTORY:
!> @date _Sep, 2018_ Initial release
!
 SUBROUTINE QGaussRys(ryscomp, cp, id, c, znuc, xyzin, igrd)
!dir$ attributes forceinline :: QGaussRys
    TYPE(shpair_t), INTENT(IN) :: cp
    INTEGER, INTENT(IN) :: id
    REAL(REAL64), INTENT(IN)   :: c(3), znuc
    REAL(REAL64), CONTIGUOUS, INTENT(OUT)  :: xyzin(0:,0:,:,:)
    INTEGER, INTENT(IN), OPTIONAL :: igrd
    type(rys_root_t), intent(inout) :: ryscomp

    INTEGER :: ni, nj, k, igrd1
    REAL(REAL64) :: ww, tt
    REAL(REAL64) :: b, d(3), dij(3)
!dir$ assume_aligned xyzin : 64

    igrd1 = 0
    IF (present(igrd)) igrd1 = igrd

    call ryscomp%evaluate()

    ASSOCIATE (pp => cp%p(id), iang => cp%iang, jang => cp%jang)
    DO k = 1, ryscomp%nroots
        ww = ryscomp%w(k)*znuc
        tt = ryscomp%u(k)/(1.0+ryscomp%u(k))
        b = 0.5*(1.0-tt)/pp%aa
        d = (pp%r-cp%rj) - tt*(pp%r-c)
        dij = cp%rj - cp%ri

        xyzin(0,0,1,k) = 1.0
        xyzin(0,0,2,k) = 1.0
        xyzin(0,0,3,k) = ww

        xyzin(1,0,1,k) = d(1)
        xyzin(1,0,2,k) = d(2)
        xyzin(1,0,3,k) = d(3)*ww

        ! VRR (Lj+1,0) <- Rpj*(Lj,0) + Lj*b*(Lj-1,0)
        DO nj = 2, (iang+jang)+igrd1
            xyzin(nj,0,:,k) = d*xyzin(nj-1,0,:,k) + (nj-1)*b*xyzin(nj-2,0,:,k)
        END DO

        ! HRR (Lj,Li+1) <- (Lj+1,Li) + Rij*(Lj,Li)
        nj = (iang+jang)+igrd1
        DO ni = 1, iang+igrd1
            nj = nj-1
            xyzin(0:nj,ni,1,k) = xyzin(1:nj+1,ni-1,1,k) + dij(1)*xyzin(0:nj,ni-1,1,k)
            xyzin(0:nj,ni,2,k) = xyzin(1:nj+1,ni-1,2,k) + dij(2)*xyzin(0:nj,ni-1,2,k)
            xyzin(0:nj,ni,3,k) = xyzin(1:nj+1,ni-1,3,k) + dij(3)*xyzin(0:nj,ni-1,3,k)
        END DO

    END DO
    END ASSOCIATE

 END SUBROUTINE

!> @brief Compute 1D Coulomb integrals with Gaussian damping
!> @details Compute 1D integrals for the modified Coulomb potential:
!>  \f$ |r-r_C|^{-1}\cdot e^{-\alpha(r-r_C)^2} \f$
!> @param[in]  nroots   roots for GaussRys
!> @param[in]  cp       shell pair data
!> @param[in]  id       current pair of primitives
!> @param[in]  c        coordinates of the charged particle
!> @param[in]  znuc     charge of the particle
!> @param[in]  alpha    dumping exponent
!> @param[out] xyzin    array of 1D integrals
!
!> @author   Vladimir Mironov
!
!     REVISION HISTORY:
!> @date _Sep, 2018_ Initial release
 SUBROUTINE QGaussRys_damp(ryscomp,cp,id,c,znuc,alpha,xyzin)
!dir$ attributes forceinline :: QGaussRys_damp
    TYPE(shpair_t), INTENT(IN) :: cp
    INTEGER, INTENT(IN) :: id
    REAL(REAL64), INTENT(IN)   :: alpha, c(3), znuc
    REAL(REAL64), CONTIGUOUS, INTENT(OUT)  :: xyzin(0:,0:,:,:)

    type(rys_root_t), intent(inout) :: ryscomp

    INTEGER :: ni, nj
    INTEGER :: k
    REAL(REAL64) :: ww, tt
    REAL(REAL64) :: b, d(3), dij(3)
!dir$ assume_aligned xyzin : 64

    call ryscomp%evaluate()

    ASSOCIATE (pp => cp%p(id), iang => cp%iang, jang => cp%jang)

    DO k = 1, ryscomp%nroots
        ww = ryscomp%w(k)*znuc
        tt = ryscomp%u(k)/(1.0+ryscomp%u(k))
!       Recurrence coefficients are slightly differend from those used in regular G-R quadrature
        b = 0.5*(1.0-tt)/(pp%aa+alpha)
        d = (c - cp%rj) + 2*pp%aa*b*(pp%r - c)
        dij = cp%rj - cp%ri

        xyzin(0,0,1,k) = 1.0
        xyzin(0,0,2,k) = 1.0
        xyzin(0,0,3,k) = ww

        xyzin(1,0,1,k) = d(1)
        xyzin(1,0,2,k) = d(2)
        xyzin(1,0,3,k) = d(3)*ww

        ! VRR (Lj+1,0) <- Rpj*(Lj,0) + Lj*b*(Lj-1,0)
        DO nj = 2, iang+jang
            xyzin(nj,0,:,k) = d(:)*xyzin(nj-1,0,:,k) + (nj-1)*b*xyzin(nj-2,0,:,k)
        END DO

        ! HRR (Lj,Li+1) <- (Lj+1,Li) + Rij*(Lj,Li)
        nj = iang+jang
        DO ni = 1, iang
            nj = nj-1
            xyzin(0:nj,ni,1,k) = xyzin(1:nj+1,ni-1,1,k) + dij(1)*xyzin(0:nj,ni-1,1,k)
            xyzin(0:nj,ni,2,k) = xyzin(1:nj+1,ni-1,2,k) + dij(2)*xyzin(0:nj,ni-1,2,k)
            xyzin(0:nj,ni,3,k) = xyzin(1:nj+1,ni-1,3,k) + dij(3)*xyzin(0:nj,ni-1,3,k)
        END DO

    END DO
    END ASSOCIATE

 END SUBROUTINE

!> @brief Compute 1D integrals for Ewald long-range 1e Coulomb integrals
!> @details In this implementation 1D integrals at Rys abscissae are
!>  computed using VRR and HRR recurrences
!> @note The common factor for the integral block is \f$ 2\Pi \f$
!> @param[in]  nroots   roots for GaussRys
!> @param[in]  cp       shell pair data
!> @param[in]  id       current pair of primitives
!> @param[in]  c        coordinates of the charged particle
!> @param[in]  znuc     charge of the particle
!> @param[in]  xfac     factor a*a/(p+a*a)
!> @param[out] xyzin    array of 1D integrals
!> @param[in]  igrd     [opt] flag indicating that integral derivatives are needed
!> @author   Vladimir Mironov
!
!     REVISION HISTORY:
!> @date _Sep, 2018_ Initial release
!
 SUBROUTINE QGaussRysEw(ryscomp, cp, id, c, znuc, xfac, xyzin, igrd)
!dir$ attributes forceinline :: QGaussRysEw
    TYPE(shpair_t), INTENT(IN) :: cp
    INTEGER, INTENT(IN) :: id
    REAL(REAL64), INTENT(IN)   :: c(3), znuc, xfac
    REAL(REAL64), CONTIGUOUS, INTENT(OUT)  :: xyzin(0:,0:,:,:)
    INTEGER, INTENT(IN), OPTIONAL :: igrd
    type(rys_root_t) :: ryscomp

    INTEGER :: ni, nj, k, igrd1
    REAL(REAL64) :: ww, tt
    REAL(REAL64) :: b, d(3), dij(3)
!dir$ assume_aligned xyzin : 64

    igrd1 = 0
    IF (present(igrd)) igrd1 = igrd

    call ryscomp%evaluate()

    ASSOCIATE (pp => cp%p(id), iang => cp%iang, jang => cp%jang)
    DO k = 1, ryscomp%nroots
        ww = ryscomp%w(k)*znuc
        tt = ryscomp%u(k)/(1.0+ryscomp%u(k))*xfac
        b = 0.5*(1.0-tt)/pp%aa
        d = (pp%r-cp%rj) - tt*(pp%r-c)
        dij = cp%rj - cp%ri

        xyzin(0,0,1,k) = 1.0
        xyzin(0,0,2,k) = 1.0
        xyzin(0,0,3,k) = ww

        xyzin(1,0,1,k) = d(1)
        xyzin(1,0,2,k) = d(2)
        xyzin(1,0,3,k) = d(3)*ww

        ! VRR (Lj+1,0) <- Rpj*(Lj,0) + Lj*b*(Lj-1,0)
        DO nj = 2, iang+jang+igrd1
            xyzin(nj,0,:,k) = d*xyzin(nj-1,0,:,k) + (nj-1)*b*xyzin(nj-2,0,:,k)
        END DO

        ! HRR (Lj,Li+1) <- (Lj+1,Li) + Rij*(Lj,Li)
        nj = iang+jang+igrd1
        DO ni = 1, iang+igrd1
            nj = nj-1
            xyzin(0:nj,ni,1,k) = xyzin(1:nj+1,ni-1,1,k) + dij(1)*xyzin(0:nj,ni-1,1,k)
            xyzin(0:nj,ni,2,k) = xyzin(1:nj+1,ni-1,2,k) + dij(2)*xyzin(0:nj,ni-1,2,k)
            xyzin(0:nj,ni,3,k) = xyzin(1:nj+1,ni-1,3,k) + dij(3)*xyzin(0:nj,ni-1,3,k)
        END DO

    END DO
    END ASSOCIATE

 END SUBROUTINE

!> @brief Compute 1D integrals needed in calculation of the  Hellmann-Feynman
!>  contribution to the gradient
!> @details They differ from the regular integrals
!>  only by the \f$ 2u^2 \f$ factor. Note, that 2*(ai+aj) factor is absent
!>  here - it will be applied to the final gradient contribution
!  TODO:
!  Redesign Gauss-Rys quadrature code to handle both cases
!> @note The common factor for the integral block is \f$ 2\Pi \f$
!> @param[in]  nroots   roots for GaussRys
!> @param[in]  cp       shell pair data
!> @param[in]  id       current pair of primitives
!> @param[in]  c        coordinates of the charged particle
!> @param[in]  znuc     charge of the particle
!> @param[out] xyzin    array of 1D integrals
!> @author   Vladimir Mironov
!
!     REVISION HISTORY:
!> @date _Sep, 2018_ Initial release
!
 SUBROUTINE DQGaussRys(ryscomp, cp, id, c, znuc, xyzin)
!dir$ attributes forceinline :: DQGaussRys
    TYPE(shpair_t), INTENT(IN) :: cp
    INTEGER, INTENT(IN) :: id
    REAL(REAL64), INTENT(IN)   :: c(3), znuc
    REAL(REAL64), CONTIGUOUS, INTENT(OUT)  :: xyzin(0:,0:,:,:)

    type(rys_root_t), intent(inout) :: ryscomp

    INTEGER :: ni, nj, k
    REAL(REAL64) :: ww, tt
    REAL(REAL64) :: b, d(3), dij(3)
!dir$ assume_aligned xyzin : 64

    call ryscomp%evaluate()

    ASSOCIATE (pp => cp%p(id), iang => cp%iang, jang => cp%jang)

    dij = cp%rj - cp%ri

    DO k = 1, ryscomp%nroots
        ww = ryscomp%w(k)*znuc*ryscomp%u(k)
        tt = ryscomp%u(k)/(1.0+ryscomp%u(k))
        b = 0.5*(1.0-tt)/pp%aa

        d = (pp%r-cp%rj) - tt*(pp%r-c)

        xyzin(0,0,1,k) = 1.0
        xyzin(0,0,2,k) = 1.0
        xyzin(0,0,3,k) = ww

        xyzin(1,0,1,k) = d(1)
        xyzin(1,0,2,k) = d(2)
        xyzin(1,0,3,k) = d(3)*ww

        ! VRR (Lj+1,0) <- Rpj*(Lj,0) + Lj*b*(Lj-1,0)
        DO nj = 2, iang+jang+1
            xyzin(nj,0,:,k) = d(:)*xyzin(nj-1,0,:,k) + (nj-1)*b*xyzin(nj-2,0,:,k)
        END DO

        ! HRR (Lj,Li+1) <- (Lj+1,Li) + Rij*(Lj,Li)
        nj = iang+jang+1
        DO ni = 1, iang+1
            nj = nj-1
            xyzin(0:nj,ni,1,k) = xyzin(1:nj+1,ni-1,1,k) + dij(1)*xyzin(0:nj,ni-1,1,k)
            xyzin(0:nj,ni,2,k) = xyzin(1:nj+1,ni-1,2,k) + dij(2)*xyzin(0:nj,ni-1,2,k)
            xyzin(0:nj,ni,3,k) = xyzin(1:nj+1,ni-1,3,k) + dij(3)*xyzin(0:nj,ni-1,3,k)
        END DO

    END DO
    END ASSOCIATE

 END SUBROUTINE

!> @brief Compute 1D integrals needed in calculation of the  Hellmann-Feynman
!>  contribution to the gradient (Ewald long-range)
!> @details They differ from the regular integrals
!>  only by the \f$ 2u^2 \f$ factor. Note, that 2*(ai+aj) factor is absent
!>  here - it will be applied to the final gradient contribution
!  TODO:
!  Redesign Gauss-Rys quadrature code to handle both cases
!> @note The common factor for the integral block is \f$ 2\Pi \f$
!> @param[in]  nroots   roots for GaussRys
!> @param[in]  cp       shell pair data
!> @param[in]  id       current pair of primitives
!> @param[in]  c        coordinates of the charged particle
!> @param[in]  znuc     charge of the particle
!> @param[in]  xfac     factor a*a/(p+a*a)
!> @param[out] xyzin    array of 1D integrals
!> @author   Vladimir Mironov
!
!     REVISION HISTORY:
!> @date _Sep, 2018_ Initial release
!
 SUBROUTINE DQGaussRysEw(ryscomp, cp, id, c, znuc, xfac, xyzin)
!dir$ attributes forceinline :: DQGaussRys
    TYPE(shpair_t), INTENT(IN) :: cp
    INTEGER, INTENT(IN) :: id
    REAL(REAL64), INTENT(IN)   :: c(3), znuc, xfac
    REAL(REAL64), CONTIGUOUS, INTENT(OUT)  :: xyzin(0:,0:,:,:)

    type(rys_root_t) :: ryscomp

    INTEGER :: ni, nj, k
    REAL(REAL64) :: ww, tt
    REAL(REAL64) :: b, d(3), dij(3)
!dir$ assume_aligned xyzin : 64

    call ryscomp%evaluate()

    ASSOCIATE (pp => cp%p(id), iang => cp%iang, jang => cp%jang)

    dij = cp%rj - cp%ri

    DO k = 1, ryscomp%nroots
        ww = ryscomp%w(k)*znuc*ryscomp%u(k)
        tt = ryscomp%u(k)/(1.0+ryscomp%u(k))*xfac
        b = 0.5*(1.0-tt)/pp%aa

        d = (pp%r-cp%rj) - tt*(pp%r-c)

        xyzin(0,0,1,k) = 1.0
        xyzin(0,0,2,k) = 1.0
        xyzin(0,0,3,k) = ww

        xyzin(1,0,1,k) = d(1)
        xyzin(1,0,2,k) = d(2)
        xyzin(1,0,3,k) = d(3)*ww

        ! VRR (Lj+1,0) <- Rpj*(Lj,0) + Lj*b*(Lj-1,0)
        DO nj = 2, iang+jang+1
            xyzin(nj,0,:,k) = d(:)*xyzin(nj-1,0,:,k) + (nj-1)*b*xyzin(nj-2,0,:,k)
        END DO

        ! HRR (Lj,Li+1) <- (Lj+1,Li) + Rij*(Lj,Li)
        nj = iang+jang+1
        DO ni = 1, iang+1
            nj = nj-1
            xyzin(0:nj,ni,1,k) = xyzin(1:nj+1,ni-1,1,k) + dij(1)*xyzin(0:nj,ni-1,1,k)
            xyzin(0:nj,ni,2,k) = xyzin(1:nj+1,ni-1,2,k) + dij(2)*xyzin(0:nj,ni-1,2,k)
            xyzin(0:nj,ni,3,k) = xyzin(1:nj+1,ni-1,3,k) + dij(3)*xyzin(0:nj,ni-1,3,k)
        END DO

    END DO
    END ASSOCIATE

 END SUBROUTINE

!--------------------------------------------------------------------------------
!       GENERAL SUPPLEMENTARY ROUTINES
!--------------------------------------------------------------------------------

!> @brief Add contribution of the 1e-integral block to the triangular matrix
!> @param[in]       shi     first shell data
!> @param[in]       shj     second shell data
!> @param[in]       mblk    square block of 1e integrals passed as 1D array
!> @param[inout]    m       packed triangular matrix of 1e integral contribution
!
!> @author   Vladimir Mironov
!
!     REVISION HISTORY:
!> @date _Sep, 2018_ Initial release
!
 SUBROUTINE update_triang_matrix(shi, shj, mblk, m)
    TYPE(shell_t), INTENT(IN) :: shi, shj
    REAL(REAL64), CONTIGUOUS, INTENT(INOUT) :: m(:)
    REAL(REAL64), CONTIGUOUS, INTENT(IN) :: mblk(:)

    INTEGER :: i, j, nn, li, lj, mj, mi, jmax
    LOGICAL :: iandj
!dir$ assume_aligned mblk : 64

    iandj = shi%shid==shj%shid

    jmax = shj%nao-1
    nn = 0
    DO i = 0, shi%nao-1
        li = shi%locao+i
        mi = (li*(li-1))/2
        IF (iandj) jmax = i
        DO j = 0, jmax
            lj = shj%locao+j
            mj = lj+mi
            nn = nn+1
            m(mj) = m(mj) + mblk(nn)
        END DO
    END DO
 END SUBROUTINE

!> @brief Add contribution of the 1e-integral block to the rectangular matrix
!> @param[in]       shi     first shell data
!> @param[in]       shj     second shell data
!> @param[in]       mblk    square block of 1e integrals passed as 1D array
!> @param[inout]    m       rectangular matrix of 1-e integral contribution
!
!> @author   Igor S. Gerasimov
!
!     REVISION HISTORY:
!> @date _Oct, 2022_ Initial release
!
 SUBROUTINE update_rectangular_matrix(shi, shj, mblk, m)
    TYPE(shell_t), INTENT(IN) :: shi, shj
    REAL(REAL64), CONTIGUOUS, INTENT(INOUT) :: m(:,:)
    REAL(REAL64), CONTIGUOUS, INTENT(IN) :: mblk(:)

    INTEGER :: i, j, nn, li, lj
!dir$ assume_aligned mblk : 64

    nn = 0
    DO i = 0, shi%nao-1
        li = shi%locao+i
        DO j = 0, shj%nao-1
            lj = shj%locao+j
            nn = nn+1
            m(lj, li) = m(lj, li) + mblk(nn)
        END DO
    END DO
 END SUBROUTINE

!> @brief Copy density block from the triangular density matrix
!> @details This subroutine assumes arbitrary order of shell IDs
!> @param[in]       shi     first shell data
!> @param[in]       shj     second shell data
!> @param[in]       dij     density matrix in packed triangular form
!> @param[out]      denab   density matrix block for shells shi and shj
!> @note Used in TVDER-based subroutines
!
!> @author   Vladimir Mironov
!
!     REVISION HISTORY:
!> @date _Sep, 2018_ Initial release
!
 SUBROUTINE density_unordered(shi, shj, dij, denab)
    TYPE(shell_t), INTENT(IN) :: shi, shj
    REAL(REAL64), CONTIGUOUS, INTENT(IN) :: denab(:)
    REAL(REAL64), CONTIGUOUS, INTENT(OUT) :: dij(:)

    INTEGER :: ij, i, j, nn, i0, j0

    ij = 0
    DO i = 0, shi%nao-1
        DO j = 0, shj%nao-1

            ij = ij+1
            i0 = max(shi%locao+i,shj%locao+j)
            j0 = min(shi%locao+i,shj%locao+j)
            nn = (i0-1)*i0/2 + j0
            dij(ij) = 2*denab(nn)
        END DO
    END DO

 END SUBROUTINE

!> @brief Copy density block from the triangular density matrix
!> @details This subroutine assumes `shi%shid>=shj%shid`
!> @param[in]       shi     first shell data
!> @param[in]       shj     second shell data
!> @param[in]       dij     density matrix in packed triangular form
!> @param[out]      denab   density matrix block for shells shi and shj
!> @note Used in HELFEY-based subroutines
!
!> @author   Vladimir Mironov
!
!     REVISION HISTORY:
!> @date _Sep, 2018_ Initial release
!
 SUBROUTINE density_ordered(shi, shj, dij, denab)
    TYPE(shell_t), INTENT(IN) :: shi, shj
    REAL(REAL64), CONTIGUOUS, INTENT(IN) :: denab(:)
    REAL(REAL64), CONTIGUOUS, INTENT(OUT) :: dij(:)

    INTEGER :: ij, i, j, nn, jmax, i0, j0
    REAL(REAL64) :: den
    LOGICAL :: iandj

    iandj = shi%shid == shj%shid

    jmax = shj%nao-1

    ij = 0
    DO i = 0, shi%nao-1

        IF (iandj) jmax = i
        DO j = 0, jmax
            ij = ij+1
            i0 = shi%locao+i
            j0 = shj%locao+j
            nn = (i0-1)*i0/2 + j0

            den = 2*denab(nn)
            IF (iandj.AND.i==j) den = denab(nn)
            dij(ij) = den

        END DO
    END DO
 END SUBROUTINE

END MODULE
