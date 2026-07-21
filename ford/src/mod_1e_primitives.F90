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
 PUBLIC comp_amom_int1_prim
 PUBLIC comp_giao_overlap_deriv_prim
 PUBLIC comp_giao_h10_core_prim
 PUBLIC comp_nmr_dia_int1_prim

 PUBLIC comp_pso_int1_prim
 PUBLIC comp_giao_a01gp_prim
 public comp_mult_int1_prim
 public comp_allmult_int1_prim
 PUBLIC comp_coulomb_dampch_int1_prim
 PUBLIC comp_ewaldlr_int1_prim
 PUBLIC comp_coulpot_prim
 PUBLIC comp_coulomb_der1
 PUBLIC comp_coulomb_helfeyder1
 PUBLIC comp_kinetic_der1
 PUBLIC comp_overlap_der1
 PUBLIC comp_overlap_der1_block
 PUBLIC comp_kinetic_der1_block
 PUBLIC comp_coulomb_der1_block
 PUBLIC comp_coulomb_helfeyder1_block
 PUBLIC comp_kinetic_der2
 PUBLIC comp_overlap_der2
 PUBLIC der_kinovl_xyz
 PUBLIC der2_kinovl_xyz
 PUBLIC der_coul_xyz
 PUBLIC der2_coul_xyz
 PUBLIC comp_coulomb_der2_braC
 PUBLIC comp_coulomb_der2_blocks
 PUBLIC comp_ewaldlr_der1
 PUBLIC comp_ewaldlr_helfeyder1

 PUBLIC update_triang_matrix
 PUBLIC update_rectangular_matrix
 PUBLIC density_ordered
 PUBLIC density_unordered
 PUBLIC comp_pvp_int1_prim
 PUBLIC comp_soc_int1_prim
 public comp_soc_int2_prim
 public QGaussRys2e

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

!> @brief Compute primitive block of angular momentum integrals about a
!>        gauge origin `o`, all three components.
!> @details Computes the real, antisymmetric matrix elements of the orbital
!>  angular momentum operator measured about the point `o`:
!>    \f$ A_a = \epsilon_{abc} (r-o)_b \partial_c \f$, a,b,c = x,y,z
!>  i.e. \f$ L = -i\,(r-o)\times\nabla = -i\,A \f$, so the physical angular
!>  momentum operator is \f$ -i \f$ times the block returned here.
!>  Each component factorizes into a product of 1D factors: a moment-about-`o`
!>  factor on one axis, a ket-derivative factor on another, and a plain overlap
!>  on the third. The ket derivative uses the Gaussian rule
!>    \f$ \partial \phi_\nu = 2\alpha_\nu \phi_{\nu+1} - n_\nu \phi_{\nu-1} \f$.
!> @param[in]       cp          shell pair data
!> @param[in]       id          current pair of primitives
!> @param[in]       o           gauge origin
!> @param[inout]    blk         block of 1e angular momentum integrals (:,1:3)
!
!> @author   Generated for NMR shielding (CGO)
!
 SUBROUTINE comp_amom_int1_prim(cp, id, o, blk)
!dir$ attributes inline :: comp_amom_int1_prim

    type(shpair_t), intent(in)  :: cp
    integer, intent(in) :: id
    real(real64), contiguous, intent(in) :: o(:)
    real(real64), contiguous, intent(inout) :: blk(:,:)

    integer, parameter :: X__ = 1, Y__ = 2, Z__ = 3
    INTEGER :: i, j, nx, ny, nz, mx, my, mz, ij, jmax
    REAL(REAL64) :: aj
    REAL(REAL64) :: sx0, sy0, sz0, mx1, my1, mz1, dx, dy, dz
    real(real64) :: xyzmom(3,0:1,0:max_ang+1,0:max_ang)

    ASSOCIATE (pp => cp%p(id))
    aj = pp%aj
    ! moments about `o` (mom 0 = overlap, mom 1 = (q-o) moment),
    ! ket angular momentum extended by 1 to allow the ket derivative.
    CALL multipole_xyz(cp%ri, cp%rj, pp%r, pp%aa1, cp%iang, cp%jang+1, o, 1, xyzmom)

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

        ! plain overlaps
        sx0 = xyzmom(X__,0,mx,nx)
        sy0 = xyzmom(Y__,0,my,ny)
        sz0 = xyzmom(Z__,0,mz,nz)
        ! (q-o) moments
        mx1 = xyzmom(X__,1,mx,nx)
        my1 = xyzmom(Y__,1,my,ny)
        mz1 = xyzmom(Z__,1,mz,nz)
        ! ket derivatives: 2*aj*S(ket+1) - n_ket*S(ket-1)
        dx = 2*aj*xyzmom(X__,0,mx+1,nx) - mx*xyzmom(X__,0,max(mx-1,0),nx)
        dy = 2*aj*xyzmom(Y__,0,my+1,ny) - my*xyzmom(Y__,0,max(my-1,0),ny)
        dz = 2*aj*xyzmom(Z__,0,mz+1,nz) - mz*xyzmom(Z__,0,max(mz-1,0),nz)

        ! A_x = (y-o_y) d_z - (z-o_z) d_y
        blk(ij,X__) = blk(ij,X__) + pp%expfac * sx0 * (my1*dz - mz1*dy)
        ! A_y = (z-o_z) d_x - (x-o_x) d_z
        blk(ij,Y__) = blk(ij,Y__) + pp%expfac * sy0 * (mz1*dx - mx1*dz)
        ! A_z = (x-o_x) d_y - (y-o_y) d_x
        blk(ij,Z__) = blk(ij,Z__) + pp%expfac * sz0 * (mx1*dy - my1*dx)
      END DO
    END DO

    END ASSOCIATE

 END SUBROUTINE

!> @brief Primitive GIAO/London overlap magnetic derivative block.
!> @details Accumulates the real coefficient of the imaginary first magnetic
!>  derivative of the AO overlap matrix, omitting the common factor i.  For a
!>  bra function centered at R_mu and ket function centered at R_nu,
!>    S10_a(mu,nu) = 0.5 * [(R_mu - R_nu) x <mu|r|nu>]_a.
!>  This is the first native GIAO one-electron building block; it is not wired
!>  into the NMR shielding dispatch until h10/two-electron/CPHF terms pass the
!>  benchmark matrix.
 SUBROUTINE comp_giao_overlap_deriv_prim(cp, id, blk)
!dir$ attributes inline :: comp_giao_overlap_deriv_prim

    type(shpair_t), intent(in)  :: cp
    integer, intent(in) :: id
    real(real64), contiguous, intent(inout) :: blk(:,:)

    integer, parameter :: X__ = 1, Y__ = 2, Z__ = 3
    INTEGER :: i, j, nx, ny, nz, mx, my, mz, ij, jmax
    REAL(REAL64) :: sx0, sy0, sz0, mx1, my1, mz1, mu_x, mu_y, mu_z
    REAL(REAL64) :: dr(3), zero(3)
    real(real64) :: xyzmom(3,0:1,0:max_ang,0:max_ang)

    zero = 0.0_real64
    dr = cp%ri(:3) - cp%rj(:3)

    ASSOCIATE (pp => cp%p(id))
    CALL multipole_xyz(cp%ri, cp%rj, pp%r, pp%aa1, cp%iang, cp%jang, zero, 1, xyzmom)

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

        sx0 = xyzmom(X__,0,mx,nx)
        sy0 = xyzmom(Y__,0,my,ny)
        sz0 = xyzmom(Z__,0,mz,nz)
        mx1 = xyzmom(X__,1,mx,nx)
        my1 = xyzmom(Y__,1,my,ny)
        mz1 = xyzmom(Z__,1,mz,nz)
        mu_x = mx1*sy0*sz0
        mu_y = sx0*my1*sz0
        mu_z = sx0*sy0*mz1

        blk(ij,X__) = blk(ij,X__) + 0.5_real64*pp%expfac*(dr(Y__)*mu_z - dr(Z__)*mu_y)
        blk(ij,Y__) = blk(ij,Y__) + 0.5_real64*pp%expfac*(dr(Z__)*mu_x - dr(X__)*mu_z)
        blk(ij,Z__) = blk(ij,Z__) + 0.5_real64*pp%expfac*(dr(X__)*mu_y - dr(Y__)*mu_x)
      END DO
    END DO

    END ASSOCIATE

 END SUBROUTINE

!> @brief Primitive GIAO/London first-order core-Hamiltonian magnetic derivative.
!> @details Accumulates the real coefficient of the imaginary RHF GIAO h10 one-
!>  electron operator, omitting the common factor i.  The convention follows the
!>  libcint RHF NMR core-orbital convention
!>  h10_core = - int1e_ignuc(asym) - int1e_igkin.  The assembled
!>  int1_giao_h10_core routine adds the separate -0.5*int1e_giao_irjxp
!>  one-electron GIAO term.  This is still not a shielding, does not include the
!>  GIAO two-electron Fock derivative, and must not ungate nmr_gauge=giao by
!>  itself.
 SUBROUTINE comp_giao_h10_core_prim(cp, id, coord, zq, nat, blk)
!dir$ attributes inline :: comp_giao_h10_core_prim

    type(shpair_t), intent(in)  :: cp
    integer, intent(in) :: id
    real(real64), contiguous, intent(in) :: coord(:,:), zq(:)
    integer, intent(in) :: nat
    real(real64), contiguous, intent(inout) :: blk(:,:)

    integer, parameter :: X__ = 1, Y__ = 2, Z__ = 3
    integer, parameter :: NRT = MAX_NROOTS
    type(rys_root_t) :: ryscomp
    integer :: i, j, ic, nx, ny, nz, mx, my, mz, ij, jmax
    real(real64) :: xx, dij, kin0, nuc0, kin_mom(3), nuc_mom(3), mom(3), cvec(3)
    real(real64) :: ovl_x0, ovl_y0, ovl_z0, ovl_x1, ovl_y1, ovl_z1
    real(real64) :: rys_x0, rys_y0, rys_z0, rys_x1, rys_y1, rys_z1
    real(real64) :: xyzovl(0:max_ang+2,0:max_ang+1,3)
    real(real64) :: xyzkin(0:max_ang_pad,0:max_ang+1,3)
    real(real64) :: xyzin(0:2*max_ang+1, 0:max_ang+1,3,NRT)
!dir$ assume_aligned blk : 64
!dir$ assume_aligned xyzkin : 64
!dir$ assume_aligned xyzovl : 64
!dir$ assume_aligned xyzin : 64

    cvec = cp%ri(:3) - cp%rj(:3)

    associate (pp => cp%p(id))
    call overlap_xyz(cp%ri, cp%rj, pp%r, pp%aa1, cp%iang+1, cp%jang+2, xyzovl)
    call kinetic_xyz_j(xyzkin, xyzovl, cp%iang+1, cp%jang, pp%aj)

    ij = 0
    jmax = cp%jnao
    do i = 1, cp%inao
      nx = CART_X(i,cp%iang)
      ny = CART_Y(i,cp%iang)
      nz = CART_Z(i,cp%iang)
      if (cp%iandj) jmax = i
      do j = 1, jmax
        mx = CART_X(j,cp%jang)
        my = CART_Y(j,cp%jang)
        mz = CART_Z(j,cp%jang)
        ij = ij + 1

        ovl_x0 = xyzovl(mx,nx,X__)
        ovl_y0 = xyzovl(my,ny,Y__)
        ovl_z0 = xyzovl(mz,nz,Z__)
        ovl_x1 = xyzovl(mx,nx+1,X__)
        ovl_y1 = xyzovl(my,ny+1,Y__)
        ovl_z1 = xyzovl(mz,nz+1,Z__)

        kin0 = xyzkin(mx,nx,X__)*ovl_y0*ovl_z0 &
             + ovl_x0*xyzkin(my,ny,Y__)*ovl_z0 &
             + ovl_x0*ovl_y0*xyzkin(mz,nz,Z__)
        kin_mom(X__) = xyzkin(mx,nx+1,X__)*ovl_y0*ovl_z0 &
             + ovl_x1*xyzkin(my,ny,Y__)*ovl_z0 &
             + ovl_x1*ovl_y0*xyzkin(mz,nz,Z__) &
             + cp%ri(X__)*kin0
        kin_mom(Y__) = xyzkin(mx,nx,X__)*ovl_y1*ovl_z0 &
             + ovl_x0*xyzkin(my,ny+1,Y__)*ovl_z0 &
             + ovl_x0*ovl_y1*xyzkin(mz,nz,Z__) &
             + cp%ri(Y__)*kin0
        kin_mom(Z__) = xyzkin(mx,nx,X__)*ovl_y0*ovl_z1 &
             + ovl_x0*xyzkin(my,ny,Y__)*ovl_z1 &
             + ovl_x0*ovl_y0*xyzkin(mz,nz+1,Z__) &
             + cp%ri(Z__)*kin0
        nuc_mom = 0.0_real64
        nuc0 = 0.0_real64

        do ic = 1, nat
          xx = pp%aa*sum((pp%r(:3) - coord(:,ic))**2)
          ryscomp%nroots = cp%nroots
          ryscomp%x = xx
          call QGaussRys(ryscomp, cp, id, coord(:,ic), -zq(ic), xyzin, 1)
          dij = pp%expfac*TWOPI*pp%aa1
          nuc0 = nuc0 + dij*sum(xyzin(mx,nx,X__,1:cp%nroots) &
                                * xyzin(my,ny,Y__,1:cp%nroots) &
                                * xyzin(mz,nz,Z__,1:cp%nroots))
          nuc_mom(X__) = nuc_mom(X__) + dij*sum(xyzin(mx,nx+1,X__,1:cp%nroots) &
                                * xyzin(my,ny,Y__,1:cp%nroots) &
                                * xyzin(mz,nz,Z__,1:cp%nroots))
          nuc_mom(Y__) = nuc_mom(Y__) + dij*sum(xyzin(mx,nx,X__,1:cp%nroots) &
                                * xyzin(my,ny+1,Y__,1:cp%nroots) &
                                * xyzin(mz,nz,Z__,1:cp%nroots))
          nuc_mom(Z__) = nuc_mom(Z__) + dij*sum(xyzin(mx,nx,X__,1:cp%nroots) &
                                * xyzin(my,ny,Y__,1:cp%nroots) &
                                * xyzin(mz,nz+1,Z__,1:cp%nroots))
        end do
        nuc_mom(:) = nuc_mom(:) + cp%ri(:)*nuc0

        mom = pp%expfac*kin_mom + nuc_mom
        blk(ij,X__) = blk(ij,X__) + 0.5_real64*(cvec(Y__)*mom(Z__) - cvec(Z__)*mom(Y__))
        blk(ij,Y__) = blk(ij,Y__) + 0.5_real64*(cvec(Z__)*mom(X__) - cvec(X__)*mom(Z__))
        blk(ij,Z__) = blk(ij,Z__) + 0.5_real64*(cvec(X__)*mom(Y__) - cvec(Y__)*mom(X__))
      end do
    end do

    end associate

 END SUBROUTINE


!> @brief Density-contracted NMR diamagnetic shielding integrals for one nucleus.
!> @details Accumulates the nine components
!>    g_ab = sum_{mu,nu} D_{mu,nu} <mu| (r-o)_a (r-c)_b / |r-c|^3 |nu>
!>  for a given nucleus at `c` and gauge origin `o`, summing over the primitive
!>  pairs of the contracted shell pair. The diamagnetic shielding tensor is then
!>    sigma^dia_{ts}(N) = (alpha^2/2) [ delta_ts * (g_xx+g_yy+g_zz) - g_{s,t} ].
!>  The field factor (r-c) and the gauge-moment factor (r-o) are both inserted by
!>  Cartesian index raising on the Rys nuclear-attraction kernel (the same
!>  mechanism as der_helfey_xyz), reusing one extra bra order for the field and a
!>  second for the (r-o) moment on the diagonal components.
!> @param[in]       cp     shell pair data
!> @param[in]       c      nucleus coordinates
!> @param[in]       o      gauge origin
!> @param[in]       den    density matrix block (i=bra, j=ket)
!> @param[inout]    gdia   3x3 accumulator for the contracted integrals
!
!> @author   Generated for NMR shielding (CGO)
!
 SUBROUTINE comp_nmr_dia_int1_prim(cp, c, o, den, gdia)
    TYPE(shpair_t), INTENT(IN) :: cp
    REAL(REAL64), INTENT(IN) :: c(3), o(3)
    REAL(REAL64), INTENT(IN) :: den(:,:)
    REAL(REAL64), INTENT(INOUT) :: gdia(3,3)

    integer, parameter :: NRT = MAX_NROOTS+3
    type(rys_root_t) :: ryscomp
    REAL(REAL64) :: xx, ww, tt, bb, dd(3), rji(3), ric(3), rio(3), fac
    INTEGER :: id, k, nr, ni, nj, a, b, m, kd(3)
    INTEGER :: i, j, ix, iy, iz, jx, jy, jz
    REAL(REAL64) :: prod, accum
    ! Rys kernel: (ket, bra, coord, root); bra extended by 2
    real(real64) :: xyzin(0:2*max_ang+3, 0:max_ang+2, 3, NRT)
    ! field factor (r-c): bra extended by 1
    real(real64) :: fld(0:max_ang, 0:max_ang+1, 3, NRT)
    ! per-coord factors by kind: 0=plain,1=field,2=moment,3=moment*field
    real(real64) :: facK(0:max_ang, 0:max_ang, 3, 0:3, NRT)
!dir$ assume_aligned xyzin : 64

    ric = cp%ri(:3) - c(:3)
    rio = cp%ri(:3) - o(:3)
    rji = cp%rj(:3) - cp%ri(:3)

    DO id = 1, cp%numpairs
      ASSOCIATE (pp => cp%p(id), &
                 iang => cp%iang, jang => cp%jang, &
                 inao => cp%inao, jnao => cp%jnao)

      nr = cp%nroots + 1
      xx = pp%aa*sum((pp%r-c)**2)
      ryscomp%nroots = nr
      ryscomp%x = xx
      call ryscomp%evaluate()

      DO k = 1, nr
        ww = ryscomp%w(k)*ryscomp%u(k)
        tt = ryscomp%u(k)/(1.0d0+ryscomp%u(k))
        bb = 0.5d0*(1.0d0-tt)/pp%aa
        dd = (pp%r-cp%rj) - tt*(pp%r-c)

        xyzin(0,0,1,k) = 1.0d0
        xyzin(0,0,2,k) = 1.0d0
        xyzin(0,0,3,k) = ww
        xyzin(1,0,1,k) = dd(1)
        xyzin(1,0,2,k) = dd(2)
        xyzin(1,0,3,k) = dd(3)*ww

        ! VRR (Lj+1,0)
        DO nj = 2, (iang+jang)+2
          xyzin(nj,0,:,k) = dd*xyzin(nj-1,0,:,k) + (nj-1)*bb*xyzin(nj-2,0,:,k)
        END DO
        ! HRR (Lj,Li+1), bra up to iang+2
        nj = (iang+jang)+2
        DO ni = 1, iang+2
          nj = nj-1
          xyzin(0:nj,ni,1,k) = xyzin(1:nj+1,ni-1,1,k) + rji(1)*xyzin(0:nj,ni-1,1,k)
          xyzin(0:nj,ni,2,k) = xyzin(1:nj+1,ni-1,2,k) + rji(2)*xyzin(0:nj,ni-1,2,k)
          xyzin(0:nj,ni,3,k) = xyzin(1:nj+1,ni-1,3,k) + rji(3)*xyzin(0:nj,ni-1,3,k)
        END DO

        ! field factor (r-c)_coord, bra 0:iang+1
        DO m = 1, 3
          fld(0:jang,0:iang+1,m,k) = xyzin(0:jang,1:iang+2,m,k) &
                                   + ric(m)*xyzin(0:jang,0:iang+1,m,k)
        END DO

        ! per-coord factors, bra 0:iang
        DO m = 1, 3
          ! kind 0: plain
          facK(0:jang,0:iang,m,0,k) = xyzin(0:jang,0:iang,m,k)
          ! kind 1: field
          facK(0:jang,0:iang,m,1,k) = fld(0:jang,0:iang,m,k)
          ! kind 2: (r-o) moment of plain
          facK(0:jang,0:iang,m,2,k) = xyzin(0:jang,1:iang+1,m,k) &
                                    + rio(m)*xyzin(0:jang,0:iang,m,k)
          ! kind 3: (r-o) moment of field  = field with extra (r-o) raise
          facK(0:jang,0:iang,m,3,k) = fld(0:jang,1:iang+1,m,k) &
                                    + rio(m)*fld(0:jang,0:iang,m,k)
        END DO
      END DO

      fac = pp%expfac*TWOPI*2.0d0

      DO a = 1, 3
        DO b = 1, 3
          ! kind per coordinate for this (a,b)
          DO m = 1, 3
            if (m==a .and. m==b) then
              kd(m) = 3
            else if (m==a) then
              kd(m) = 2
            else if (m==b) then
              kd(m) = 1
            else
              kd(m) = 0
            end if
          END DO

          accum = 0.0d0
          DO i = 1, inao
            ix = CART_X(i,iang); iy = CART_Y(i,iang); iz = CART_Z(i,iang)
            DO j = 1, jnao
              jx = CART_X(j,jang); jy = CART_Y(j,jang); jz = CART_Z(j,jang)
              prod = sum( facK(jx,ix,1,kd(1),1:nr) &
                        * facK(jy,iy,2,kd(2),1:nr) &
                        * facK(jz,iz,3,kd(3),1:nr) )
              accum = accum + den(i,j)*prod
            END DO
          END DO
          gdia(a,b) = gdia(a,b) + fac*accum
        END DO
      END DO

      END ASSOCIATE
    END DO

 END SUBROUTINE

!> @brief Compute primitive block of PSO (paramagnetic spin-orbit) integrals
!>        for one nucleus, all three components.
!> @details Real antisymmetric matrix elements of the operator
!>    A_a = [(r-c) x grad]_a / |r-c|^3     (so the physical PSO operator is -i*A).
!>  The field factor (r-c)/|r-c|^3 is inserted by a der_helfey-style raise on the
!>  Rys nuclear kernel (validated in the diamagnetic term); the grad factor is the
!>  ket derivative 2*aj*S(ket+1) - n_ket*S(ket-1) (as in the angular momentum term).
!>  Loops over the primitive pairs of the shell pair internally.
!> @param[in]       cp     shell pair data
!> @param[in]       c      nucleus coordinates
!> @param[inout]    blk    block of PSO integrals (:,1:3)
!
!> @author   Generated for NMR shielding (CGO)
!
 SUBROUTINE comp_pso_int1_prim(cp, c, blk)
    TYPE(shpair_t), INTENT(IN) :: cp
    REAL(REAL64), INTENT(IN) :: c(3)
    REAL(REAL64), CONTIGUOUS, INTENT(INOUT) :: blk(:,:)

    integer, parameter :: X__ = 1, Y__ = 2, Z__ = 3
    integer, parameter :: NRT = MAX_NROOTS+3
    type(rys_root_t) :: ryscomp
    REAL(REAL64) :: xx, ww, tt, bb, dd(3), rji(3), ric(3), fac, aj
    INTEGER :: id, k, nr, ni, nj, m
    INTEGER :: i, j, ix, iy, iz, jx, jy, jz, ij, jmax
    REAL(REAL64) :: px, py, pz
    real(real64) :: xyzin(0:2*max_ang+3, 0:max_ang+2, 3, NRT)
    real(real64) :: fld(0:max_ang, 0:max_ang, 3, NRT)   ! (ket,bra,coord,root)
    real(real64) :: dkt(0:max_ang, 0:max_ang, 3, NRT)   ! ket derivative
!dir$ assume_aligned xyzin : 64

    ric = cp%ri(:3) - c(:3)
    rji = cp%rj(:3) - cp%ri(:3)

    DO id = 1, cp%numpairs
      ASSOCIATE (pp => cp%p(id), &
                 iang => cp%iang, jang => cp%jang, &
                 inao => cp%inao, jnao => cp%jnao)

      aj = pp%aj
      nr = cp%nroots + 1
      xx = pp%aa*sum((pp%r-c)**2)
      ryscomp%nroots = nr
      ryscomp%x = xx
      call ryscomp%evaluate()

      DO k = 1, nr
        ww = ryscomp%w(k)*ryscomp%u(k)
        tt = ryscomp%u(k)/(1.0d0+ryscomp%u(k))
        bb = 0.5d0*(1.0d0-tt)/pp%aa
        dd = (pp%r-cp%rj) - tt*(pp%r-c)

        xyzin(0,0,1,k) = 1.0d0
        xyzin(0,0,2,k) = 1.0d0
        xyzin(0,0,3,k) = ww
        xyzin(1,0,1,k) = dd(1)
        xyzin(1,0,2,k) = dd(2)
        xyzin(1,0,3,k) = dd(3)*ww

        DO nj = 2, (iang+jang)+2
          xyzin(nj,0,:,k) = dd*xyzin(nj-1,0,:,k) + (nj-1)*bb*xyzin(nj-2,0,:,k)
        END DO
        nj = (iang+jang)+2
        DO ni = 1, iang+1
          nj = nj-1
          xyzin(0:nj,ni,1,k) = xyzin(1:nj+1,ni-1,1,k) + rji(1)*xyzin(0:nj,ni-1,1,k)
          xyzin(0:nj,ni,2,k) = xyzin(1:nj+1,ni-1,2,k) + rji(2)*xyzin(0:nj,ni-1,2,k)
          xyzin(0:nj,ni,3,k) = xyzin(1:nj+1,ni-1,3,k) + rji(3)*xyzin(0:nj,ni-1,3,k)
        END DO

        ! field factor (r-c)_coord, bra 0:iang
        DO m = 1, 3
          fld(0:jang,0:iang,m,k) = xyzin(0:jang,1:iang+1,m,k) &
                                 + ric(m)*xyzin(0:jang,0:iang,m,k)
        END DO
        ! ket derivative (2*aj*S(ket+1) - n_ket*S(ket-1)), ket 0:jang, bra 0:iang
        dkt(0,0:iang,1:3,k) = 2.0d0*aj*xyzin(1,0:iang,1:3,k)
        DO nj = 1, jang
          dkt(nj,0:iang,1:3,k) = 2.0d0*aj*xyzin(nj+1,0:iang,1:3,k) &
                               - nj*xyzin(nj-1,0:iang,1:3,k)
        END DO
      END DO

      ! The raw field+ket-derivative product carries a small spurious symmetric
      ! component for off-center nuclei. The full (non-packed) block is emitted
      ! here so the caller (pso_integrals) can antisymmetrise A=(M-M^T)/2, which
      ! is exact for the anti-Hermitian PSO operator and removes that error.
      ! Hence NO iandj triangular packing below.
      fac = pp%expfac*TWOPI*2.0d0

      ij = 0
      jmax = jnao
      DO i = 1, inao
        ix = CART_X(i,iang); iy = CART_Y(i,iang); iz = CART_Z(i,iang)
        DO j = 1, jmax
          jx = CART_X(j,jang); jy = CART_Y(j,jang); jz = CART_Z(j,jang)
          ij = ij+1

          ! PSO_x = P_x (fld_y dz - dy fld_z); etc (P=plain xyzin)
          px = sum( xyzin(jx,ix,1,1:nr) * &
                    ( fld(jy,iy,2,1:nr)*dkt(jz,iz,3,1:nr) &
                    - dkt(jy,iy,2,1:nr)*fld(jz,iz,3,1:nr) ) )
          py = sum( xyzin(jy,iy,2,1:nr) * &
                    ( fld(jz,iz,3,1:nr)*dkt(jx,ix,1,1:nr) &
                    - dkt(jz,iz,3,1:nr)*fld(jx,ix,1,1:nr) ) )
          pz = sum( xyzin(jz,iz,3,1:nr) * &
                    ( fld(jx,ix,1,1:nr)*dkt(jy,iy,2,1:nr) &
                    - dkt(jx,ix,1,1:nr)*fld(jy,iy,2,1:nr) ) )

          blk(ij,X__) = blk(ij,X__) + fac*px
          blk(ij,Y__) = blk(ij,Y__) + fac*py
          blk(ij,Z__) = blk(ij,Z__) + fac*pz
        END DO
      END DO

      END ASSOCIATE
    END DO

 END SUBROUTINE

!> @brief Primitive GIAO a01gp gauge-correction integrals (9 components).
!> @details a01gp = (g | nabla-rinv cross p |): the GIAO/London first-order
!>  derivative of the PSO operator at nucleus c.  Returns
!>    blk(ij,(a-1)*3+col) = (cvec x M^{(col)})_a ,
!>  with cvec = R_bra - R_ket and M^{(col)}_b = <mu|(r-R_bra)_b PSO_col|nu>
!>  (the bra-position-weighted PSO, built by raising the bra angular momentum by
!>  one in coordinate b, no center shift).  Full (both-triangle) block, no
!>  packing.  The overall sign/scale is calibrated by the caller against the
!>  libcint int1e_a01gp oracle.
 SUBROUTINE comp_giao_a01gp_prim(cp, c, cvec, blk)
    TYPE(shpair_t), INTENT(IN) :: cp
    REAL(REAL64), INTENT(IN) :: c(3), cvec(3)
    REAL(REAL64), CONTIGUOUS, INTENT(INOUT) :: blk(:,:)

    integer, parameter :: X__ = 1, Y__ = 2, Z__ = 3
    integer, parameter :: NRT = MAX_NROOTS+3
    type(rys_root_t) :: ryscomp
    REAL(REAL64) :: xx, ww, tt, bb, dd(3), rji(3), ric(3), fac, aj
    INTEGER :: id, k, nr, ni, nj, m, a, col
    INTEGER :: i, j, ix, iy, iz, jx, jy, jz, ij, jmax
    REAL(REAL64) :: mm(3,3), pbase(3)   ! M^{(col)}_b : (b, col); base PSO
    real(real64) :: xyzin(0:2*max_ang+3, 0:max_ang+2, 3, NRT)
    real(real64) :: fld(0:max_ang, 0:max_ang+1, 3, NRT)   ! (ket,bra,coord,root)
    real(real64) :: dkt(0:max_ang, 0:max_ang+1, 3, NRT)   ! ket derivative
!dir$ assume_aligned xyzin : 64

    ric = cp%ri(:3) - c(:3)
    rji = cp%rj(:3) - cp%ri(:3)

    DO id = 1, cp%numpairs
      ASSOCIATE (pp => cp%p(id), &
                 iang => cp%iang, jang => cp%jang, &
                 inao => cp%inao, jnao => cp%jnao)

      aj = pp%aj
      ! One more Rys root than the PSO term: the extra bra-position raise
      ! (r-R_bra) increases the polynomial order by one.
      nr = cp%nroots + 2
      xx = pp%aa*sum((pp%r-c)**2)
      ryscomp%nroots = nr
      ryscomp%x = xx
      call ryscomp%evaluate()

      DO k = 1, nr
        ww = ryscomp%w(k)*ryscomp%u(k)
        tt = ryscomp%u(k)/(1.0d0+ryscomp%u(k))
        bb = 0.5d0*(1.0d0-tt)/pp%aa
        dd = (pp%r-cp%rj) - tt*(pp%r-c)

        xyzin(0,0,1,k) = 1.0d0
        xyzin(0,0,2,k) = 1.0d0
        xyzin(0,0,3,k) = ww
        xyzin(1,0,1,k) = dd(1)
        xyzin(1,0,2,k) = dd(2)
        xyzin(1,0,3,k) = dd(3)*ww

        DO nj = 2, (iang+jang)+3
          xyzin(nj,0,:,k) = dd*xyzin(nj-1,0,:,k) + (nj-1)*bb*xyzin(nj-2,0,:,k)
        END DO
        nj = (iang+jang)+3
        DO ni = 1, iang+2
          nj = nj-1
          xyzin(0:nj,ni,1,k) = xyzin(1:nj+1,ni-1,1,k) + rji(1)*xyzin(0:nj,ni-1,1,k)
          xyzin(0:nj,ni,2,k) = xyzin(1:nj+1,ni-1,2,k) + rji(2)*xyzin(0:nj,ni-1,2,k)
          xyzin(0:nj,ni,3,k) = xyzin(1:nj+1,ni-1,3,k) + rji(3)*xyzin(0:nj,ni-1,3,k)
        END DO

        ! field factor (r-c)_coord, bra 0:iang+1
        DO m = 1, 3
          fld(0:jang,0:iang+1,m,k) = xyzin(0:jang,1:iang+2,m,k) &
                                   + ric(m)*xyzin(0:jang,0:iang+1,m,k)
        END DO
        ! ket derivative, bra 0:iang+1
        dkt(0,0:iang+1,1:3,k) = 2.0d0*aj*xyzin(1,0:iang+1,1:3,k)
        DO nj = 1, jang
          dkt(nj,0:iang+1,1:3,k) = 2.0d0*aj*xyzin(nj+1,0:iang+1,1:3,k) &
                                 - nj*xyzin(nj-1,0:iang+1,1:3,k)
        END DO
      END DO

      fac = pp%expfac*TWOPI*2.0d0

      ij = 0
      jmax = jnao
      DO i = 1, inao
        ix = CART_X(i,iang); iy = CART_Y(i,iang); iz = CART_Z(i,iang)
        DO j = 1, jmax
          jx = CART_X(j,jang); jy = CART_Y(j,jang); jz = CART_Z(j,jang)
          ij = ij+1

          ! M^{(col)}_b : bra-position-weighted PSO (raise bra coord b by one).
          ! col=1 (PSO_x):
          mm(1,1) = sum( xyzin(jx,ix+1,1,1:nr) * &
                    ( fld(jy,iy,2,1:nr)*dkt(jz,iz,3,1:nr) - dkt(jy,iy,2,1:nr)*fld(jz,iz,3,1:nr) ) )
          mm(2,1) = sum( xyzin(jx,ix,1,1:nr) * &
                    ( fld(jy,iy+1,2,1:nr)*dkt(jz,iz,3,1:nr) - dkt(jy,iy+1,2,1:nr)*fld(jz,iz,3,1:nr) ) )
          mm(3,1) = sum( xyzin(jx,ix,1,1:nr) * &
                    ( fld(jy,iy,2,1:nr)*dkt(jz,iz+1,3,1:nr) - dkt(jy,iy,2,1:nr)*fld(jz,iz+1,3,1:nr) ) )
          ! col=2 (PSO_y):
          mm(1,2) = sum( xyzin(jy,iy,2,1:nr) * &
                    ( fld(jz,iz,3,1:nr)*dkt(jx,ix+1,1,1:nr) - dkt(jz,iz,3,1:nr)*fld(jx,ix+1,1,1:nr) ) )
          mm(2,2) = sum( xyzin(jy,iy+1,2,1:nr) * &
                    ( fld(jz,iz,3,1:nr)*dkt(jx,ix,1,1:nr) - dkt(jz,iz,3,1:nr)*fld(jx,ix,1,1:nr) ) )
          mm(3,2) = sum( xyzin(jy,iy,2,1:nr) * &
                    ( fld(jz,iz+1,3,1:nr)*dkt(jx,ix,1,1:nr) - dkt(jz,iz+1,3,1:nr)*fld(jx,ix,1,1:nr) ) )
          ! col=3 (PSO_z):
          mm(1,3) = sum( xyzin(jz,iz,3,1:nr) * &
                    ( fld(jx,ix+1,1,1:nr)*dkt(jy,iy,2,1:nr) - dkt(jx,ix+1,1,1:nr)*fld(jy,iy,2,1:nr) ) )
          mm(2,3) = sum( xyzin(jz,iz,3,1:nr) * &
                    ( fld(jx,ix,1,1:nr)*dkt(jy,iy+1,2,1:nr) - dkt(jx,ix,1,1:nr)*fld(jy,iy+1,2,1:nr) ) )
          mm(3,3) = sum( xyzin(jz,iz+1,3,1:nr) * &
                    ( fld(jx,ix,1,1:nr)*dkt(jy,iy,2,1:nr) - dkt(jx,ix,1,1:nr)*fld(jy,iy,2,1:nr) ) )

          ! Base (un-raised) PSO components, to complete R0I = (r-R_bra) + R_bra,
          ! i.e. the position is referenced to the molecular origin (libcint
          ! convention), not the bra center.  M^{(col)}_b += R_bra,b * PSO_col.
          pbase(1) = sum( xyzin(jx,ix,1,1:nr) * &
                    ( fld(jy,iy,2,1:nr)*dkt(jz,iz,3,1:nr) - dkt(jy,iy,2,1:nr)*fld(jz,iz,3,1:nr) ) )
          pbase(2) = sum( xyzin(jy,iy,2,1:nr) * &
                    ( fld(jz,iz,3,1:nr)*dkt(jx,ix,1,1:nr) - dkt(jz,iz,3,1:nr)*fld(jx,ix,1,1:nr) ) )
          pbase(3) = sum( xyzin(jz,iz,3,1:nr) * &
                    ( fld(jx,ix,1,1:nr)*dkt(jy,iy,2,1:nr) - dkt(jx,ix,1,1:nr)*fld(jy,iy,2,1:nr) ) )
          do col = 1, 3
            mm(1,col) = mm(1,col) + cp%ri(1)*pbase(col)
            mm(2,col) = mm(2,col) + cp%ri(2)*pbase(col)
            mm(3,col) = mm(3,col) + cp%ri(3)*pbase(col)
          end do

          DO col = 1, 3
            blk(ij,(1-1)*3+col) = blk(ij,(1-1)*3+col) + &
              fac*( cvec(Y__)*mm(Z__,col) - cvec(Z__)*mm(Y__,col) )
            blk(ij,(2-1)*3+col) = blk(ij,(2-1)*3+col) + &
              fac*( cvec(Z__)*mm(X__,col) - cvec(X__)*mm(Z__,col) )
            blk(ij,(3-1)*3+col) = blk(ij,(3-1)*3+col) + &
              fac*( cvec(X__)*mm(Y__,col) - cvec(Y__)*mm(X__,col) )
          END DO
        END DO
      END DO

      END ASSOCIATE
    END DO

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

!> @brief Bra-center first derivative of the overlap integral, returned as an
!>   (inao, jnao, 3) block (NOT contracted with a density). Used to assemble the
!>   AO derivative-overlap matrix dS/dR for the CPHF right-hand side.
!> @param[in]    cp     shell pair data
!> @param[inout] dblk   (inao, jnao, 3) accumulated bra-center derivatives
 SUBROUTINE comp_overlap_der1_block(cp, dblk)
    TYPE(shpair_t), INTENT(IN) :: cp
    REAL(REAL64), CONTIGUOUS, INTENT(INOUT) :: dblk(:,:,:)

    INTEGER :: i, j, k, ix, iy, iz, jx, jy, jz
    real(real64) :: ovl_int(0:max_ang,0:max_ang+3,3)
    real(real64) :: ovl_der(0:max_ang,0:max_ang,3)

    DO k = 1, cp%numpairs
    ASSOCIATE (pp => cp%p(k), &
               iang => cp%iang, jang => cp%jang, &
               inao => cp%inao, jnao => cp%jnao)
    CALL overlap_xyz(cp%ri, cp%rj, pp%r, pp%aa1, iang+1, jang, ovl_int)
    CALL der_kinovl_xyz(ovl_der, ovl_int, iang, jang, pp%ai)
    DO i = 1, inao
        ix = CART_X(i,iang); iy = CART_Y(i,iang); iz = CART_Z(i,iang)
        DO j = 1, jnao
            jx = CART_X(j,jang); jy = CART_Y(j,jang); jz = CART_Z(j,jang)
            dblk(i,j,1) = dblk(i,j,1) + ovl_der(jx,ix,1)*ovl_int(jy,iy,2)*ovl_int(jz,iz,3)*pp%expfac
            dblk(i,j,2) = dblk(i,j,2) + ovl_int(jx,ix,1)*ovl_der(jy,iy,2)*ovl_int(jz,iz,3)*pp%expfac
            dblk(i,j,3) = dblk(i,j,3) + ovl_int(jx,ix,1)*ovl_int(jy,iy,2)*ovl_der(jz,iz,3)*pp%expfac
        END DO
    END DO
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

!> @brief Bra-center first derivative of the kinetic-energy integral, returned as
!>   an (inao, jnao, 3) block (not contracted), for the dT/dR matrix used in the
!>   CPHF right-hand side. Mirrors comp_kinetic_der1.
 SUBROUTINE comp_kinetic_der1_block(cp, dblk)
    TYPE(shpair_t), INTENT(IN) :: cp
    REAL(REAL64), CONTIGUOUS, INTENT(INOUT) :: dblk(:,:,:)

    INTEGER :: i, j, k, ix, iy, iz, jx, jy, jz
    real(real64) :: ovl_int(0:max_ang,0:max_ang+3,3)
    real(real64) :: ovl_der(0:max_ang,0:max_ang,3)
    real(real64) :: kin_int(0:max_ang,0:max_ang+1,3)
    real(real64) :: kin_der(0:max_ang,0:max_ang,3)

    DO k = 1, cp%numpairs
    ASSOCIATE (pp => cp%p(k), &
               iang => cp%iang, jang => cp%jang, &
               inao => cp%inao, jnao => cp%jnao)
    CALL overlap_xyz(cp%ri, cp%rj, pp%r, pp%aa1, iang+3, jang, ovl_int)
    CALL kinetic_xyz_i(kin_int, ovl_int, iang+1, jang, pp%ai)
    CALL der_kinovl_xyz(ovl_der, ovl_int, iang, jang, pp%ai)
    CALL der_kinovl_xyz(kin_der, kin_int, iang, jang, pp%ai)
    DO i = 1, inao
        ix = CART_X(i,iang); iy = CART_Y(i,iang); iz = CART_Z(i,iang)
        DO j = 1, jnao
            jx = CART_X(j,jang); jy = CART_Y(j,jang); jz = CART_Z(j,jang)
            dblk(i,j,1) = dblk(i,j,1) + ( &
                kin_der(jx,ix,1)*ovl_int(jy,iy,2)*ovl_int(jz,iz,3) + &
                ovl_der(jx,ix,1)*kin_int(jy,iy,2)*ovl_int(jz,iz,3) + &
                ovl_der(jx,ix,1)*ovl_int(jy,iy,2)*kin_int(jz,iz,3) )*pp%expfac
            dblk(i,j,2) = dblk(i,j,2) + ( &
                kin_int(jx,ix,1)*ovl_der(jy,iy,2)*ovl_int(jz,iz,3) + &
                ovl_int(jx,ix,1)*kin_der(jy,iy,2)*ovl_int(jz,iz,3) + &
                ovl_int(jx,ix,1)*ovl_der(jy,iy,2)*kin_int(jz,iz,3) )*pp%expfac
            dblk(i,j,3) = dblk(i,j,3) + ( &
                kin_int(jx,ix,1)*ovl_int(jy,iy,2)*ovl_der(jz,iz,3) + &
                ovl_int(jx,ix,1)*kin_int(jy,iy,2)*ovl_der(jz,iz,3) + &
                ovl_int(jx,ix,1)*ovl_int(jy,iy,2)*kin_der(jz,iz,3) )*pp%expfac
        END DO
    END DO
    END ASSOCIATE
    END DO
 END SUBROUTINE

!> @brief Bra-center second derivative of the 1e overlap contribution.
!> @details Accumulates the symmetric 3x3 block of second derivatives of
!>  sum_ij dij*S_ij with respect to the bra (center i) Cartesian coordinates,
!>  i.e. d2/dA_a dA_b. The full Hessian's other blocks (A-B, B-B) follow from
!>  translational invariance of the two-center integral.
!> @param[in]    cp    shell pair data
!> @param[in]    dij   density (or energy-weighted density) matrix block
!> @param[inout] de2   dimension(3,3), accumulated bra-center 2nd derivatives
 SUBROUTINE comp_overlap_der2(cp, dij, de2)
!dir$ attributes inline :: comp_overlap_der2
    TYPE(shpair_t), INTENT(IN) :: cp
    REAL(REAL64), INTENT(IN) :: dij(:,:)
    REAL(REAL64), CONTIGUOUS, INTENT(INOUT) :: de2(:,:)

    REAL(REAL64) :: de_loc(3,3), w
    REAL(REAL64) :: sx, sy, sz, dx, dy, dz, d2x, d2y, d2z
    INTEGER :: i, j, k, ix, iy, iz, jx, jy, jz
    real(real64) :: ovl_int(0:max_ang,0:max_ang+4,3)
    real(real64) :: ovl_der(0:max_ang,0:max_ang,3)
    real(real64) :: ovl_der2(0:max_ang,0:max_ang,3)

    DO k = 1, cp%numpairs
    ASSOCIATE (pp => cp%p(k), &
               iang => cp%iang, jang => cp%jang, &
               inao => cp%inao, jnao => cp%jnao)
    ! compute overlap [i+2|j]
    CALL overlap_xyz(cp%ri, cp%rj, pp%r, pp%aa1, iang+2, jang, ovl_int)

    ! compute 1D overlap 1st and 2nd derivatives [i|j]
    CALL der_kinovl_xyz(ovl_der, ovl_int, iang, jang, pp%ai)
    CALL der2_kinovl_xyz(ovl_der2, ovl_int, iang, jang, pp%ai)

    de_loc = 0.0
    DO i = 1, inao
        ix = CART_X(i,iang)
        iy = CART_Y(i,iang)
        iz = CART_Z(i,iang)
        DO j = 1, jnao
            jx = CART_X(j,jang)
            jy = CART_Y(j,jang)
            jz = CART_Z(j,jang)
            sx = ovl_int(jx,ix,1); sy = ovl_int(jy,iy,2); sz = ovl_int(jz,iz,3)
            dx = ovl_der(jx,ix,1); dy = ovl_der(jy,iy,2); dz = ovl_der(jz,iz,3)
            d2x = ovl_der2(jx,ix,1); d2y = ovl_der2(jy,iy,2); d2z = ovl_der2(jz,iz,3)
            w = dij(i,j)
            de_loc(1,1) = de_loc(1,1) + w * d2x*sy*sz
            de_loc(2,2) = de_loc(2,2) + w * sx*d2y*sz
            de_loc(3,3) = de_loc(3,3) + w * sx*sy*d2z
            de_loc(2,1) = de_loc(2,1) + w * dx*dy*sz
            de_loc(3,1) = de_loc(3,1) + w * dx*sy*dz
            de_loc(3,2) = de_loc(3,2) + w * sx*dy*dz
        END DO
    END DO
    de2(1,1) = de2(1,1) + de_loc(1,1)*pp%expfac
    de2(2,2) = de2(2,2) + de_loc(2,2)*pp%expfac
    de2(3,3) = de2(3,3) + de_loc(3,3)*pp%expfac
    de2(2,1) = de2(2,1) + de_loc(2,1)*pp%expfac
    de2(1,2) = de2(1,2) + de_loc(2,1)*pp%expfac
    de2(3,1) = de2(3,1) + de_loc(3,1)*pp%expfac
    de2(1,3) = de2(1,3) + de_loc(3,1)*pp%expfac
    de2(3,2) = de2(3,2) + de_loc(3,2)*pp%expfac
    de2(2,3) = de2(2,3) + de_loc(3,2)*pp%expfac
    END ASSOCIATE
    END DO
 END SUBROUTINE

!> @brief Bra-center second derivative of the 1e kinetic-energy contribution.
!> @details Accumulates the symmetric 3x3 block d2/dA_a dA_b of
!>  sum_ij dij*T_ij. As with comp_kinetic_der1, the kinetic operator factorizes
!>  across the three Cartesian directions as T = Tx*Sy*Sz + Sx*Ty*Sz + Sx*Sy*Tz.
!> @param[in]    cp    shell pair data
!> @param[in]    dij   density (or energy-weighted density) matrix block
!> @param[inout] de2   dimension(3,3), accumulated bra-center 2nd derivatives
 SUBROUTINE comp_kinetic_der2(cp, dij, de2)
!dir$ attributes inline :: comp_kinetic_der2
    TYPE(shpair_t), INTENT(IN) :: cp
    REAL(REAL64), INTENT(IN) :: dij(:,:)
    REAL(REAL64), CONTIGUOUS, INTENT(INOUT) :: de2(:,:)

    REAL(REAL64) :: de_loc(3,3), w
    REAL(REAL64) :: sx, sy, sz, dsx, dsy, dsz, d2sx, d2sy, d2sz
    REAL(REAL64) :: tx, ty, tz, dtx, dty, dtz, d2tx, d2ty, d2tz
    INTEGER :: i, j, k, ix, iy, iz, jx, jy, jz
    real(real64) :: ovl_int(0:max_ang,0:max_ang+4,3)
    real(real64) :: ovl_der(0:max_ang,0:max_ang,3)
    real(real64) :: ovl_der2(0:max_ang,0:max_ang,3)
    real(real64) :: kin_int(0:max_ang,0:max_ang+2,3)
    real(real64) :: kin_der(0:max_ang,0:max_ang,3)
    real(real64) :: kin_der2(0:max_ang,0:max_ang,3)

    DO k = 1, cp%numpairs
    ASSOCIATE (pp => cp%p(k), &
               iang => cp%iang, jang => cp%jang, &
               inao => cp%inao, jnao => cp%jnao)
    ! compute overlap [i+4|j]
    CALL overlap_xyz(cp%ri, cp%rj, pp%r, pp%aa1, iang+4, jang, ovl_int)

    ! compute 1D kinetic [i+2|j]
    CALL kinetic_xyz_i(kin_int, ovl_int, iang+2, jang, pp%ai)

    ! 1st and 2nd bra-center derivatives of 1D overlap and kinetic
    CALL der_kinovl_xyz(ovl_der,  ovl_int, iang, jang, pp%ai)
    CALL der2_kinovl_xyz(ovl_der2, ovl_int, iang, jang, pp%ai)
    CALL der_kinovl_xyz(kin_der,  kin_int, iang, jang, pp%ai)
    CALL der2_kinovl_xyz(kin_der2, kin_int, iang, jang, pp%ai)

    de_loc = 0.0
    DO i = 1, inao
        ix = CART_X(i,iang)
        iy = CART_Y(i,iang)
        iz = CART_Z(i,iang)
        DO j = 1, jnao
            jx = CART_X(j,jang)
            jy = CART_Y(j,jang)
            jz = CART_Z(j,jang)
            sx = ovl_int(jx,ix,1);  sy = ovl_int(jy,iy,2);  sz = ovl_int(jz,iz,3)
            dsx = ovl_der(jx,ix,1); dsy = ovl_der(jy,iy,2); dsz = ovl_der(jz,iz,3)
            d2sx = ovl_der2(jx,ix,1); d2sy = ovl_der2(jy,iy,2); d2sz = ovl_der2(jz,iz,3)
            tx = kin_int(jx,ix,1);  ty = kin_int(jy,iy,2);  tz = kin_int(jz,iz,3)
            dtx = kin_der(jx,ix,1); dty = kin_der(jy,iy,2); dtz = kin_der(jz,iz,3)
            d2tx = kin_der2(jx,ix,1); d2ty = kin_der2(jy,iy,2); d2tz = kin_der2(jz,iz,3)
            w = dij(i,j)
            ! diagonal: d2/dA_a^2 of (Tx Sy Sz + Sx Ty Sz + Sx Sy Tz)
            de_loc(1,1) = de_loc(1,1) + w*( d2tx*sy*sz + d2sx*ty*sz + d2sx*sy*tz )
            de_loc(2,2) = de_loc(2,2) + w*( tx*d2sy*sz + sx*d2ty*sz + sx*d2sy*tz )
            de_loc(3,3) = de_loc(3,3) + w*( tx*sy*d2sz + sx*ty*d2sz + sx*sy*d2tz )
            ! off-diagonal: d2/dA_a dA_b
            de_loc(2,1) = de_loc(2,1) + w*( dtx*dsy*sz + dsx*dty*sz + dsx*dsy*tz )
            de_loc(3,1) = de_loc(3,1) + w*( dtx*sy*dsz + dsx*ty*dsz + dsx*sy*dtz )
            de_loc(3,2) = de_loc(3,2) + w*( tx*dsy*dsz + sx*dty*dsz + sx*dsy*dtz )
        END DO
    END DO
    de2(1,1) = de2(1,1) + de_loc(1,1)*pp%expfac
    de2(2,2) = de2(2,2) + de_loc(2,2)*pp%expfac
    de2(3,3) = de2(3,3) + de_loc(3,3)*pp%expfac
    de2(2,1) = de2(2,1) + de_loc(2,1)*pp%expfac
    de2(1,2) = de2(1,2) + de_loc(2,1)*pp%expfac
    de2(3,1) = de2(3,1) + de_loc(3,1)*pp%expfac
    de2(1,3) = de2(1,3) + de_loc(3,1)*pp%expfac
    de2(3,2) = de2(3,2) + de_loc(3,2)*pp%expfac
    de2(2,3) = de2(2,3) + de_loc(3,2)*pp%expfac
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

!> @brief Bra-center first derivative of the nuclear-attraction integral for one
!>   charge, returned as an (inao, jnao, 3) block (not contracted). Mirrors
!>   comp_coulomb_der1; used to assemble the dV/dR matrix for the CPHF RHS.
 SUBROUTINE comp_coulomb_der1_block(cp, c, znuc, dblk)
    TYPE(shpair_t), INTENT(IN) :: cp
    REAL(REAL64), INTENT(IN) :: c(3), znuc
    REAL(REAL64), CONTIGUOUS, INTENT(INOUT) :: dblk(:,:,:)

    REAL(REAL64) :: xx, fac, der(3)
    type(rys_root_t) :: ryscomp
    INTEGER :: id, i, j, ix, iy, iz, jx, jy, jz
    real(real64) :: xyzin(0:2*max_ang+1, 0:max_ang+1,3,max_nroots)
    real(real64) :: dxyzc(0:max_ang_pad,0:max_ang,3,max_nroots)

    DO id = 1, cp%numpairs
        ASSOCIATE (pp => cp%p(id), iang => cp%iang, jang => cp%jang, &
                   inao => cp%inao, jnao => cp%jnao)
        xx = pp%aa*sum((pp%r-c)**2)
        ryscomp%nroots = cp%nroots
        ryscomp%x = xx
        CALL QGaussRys(ryscomp, cp, id, c, znuc, xyzin, 1)
        CALL der_coul_xyz(dxyzc, xyzin, iang, jang, pp%ai, cp%nroots)
        fac = pp%expfac*TWOPI*pp%aa1
        DO i = 1, inao
            ix = CART_X(i,iang); iy = CART_Y(i,iang); iz = CART_Z(i,iang)
            DO j = 1, jnao
                jx = CART_X(j,jang); jy = CART_Y(j,jang); jz = CART_Z(j,jang)
                der(1) = sum(dxyzc(jx,ix,1,1:cp%nroots)*xyzin(jy,iy,2,1:cp%nroots)*xyzin(jz,iz,3,1:cp%nroots))
                der(2) = sum(xyzin(jx,ix,1,1:cp%nroots)*dxyzc(jy,iy,2,1:cp%nroots)*xyzin(jz,iz,3,1:cp%nroots))
                der(3) = sum(xyzin(jx,ix,1,1:cp%nroots)*xyzin(jy,iy,2,1:cp%nroots)*dxyzc(jz,iz,3,1:cp%nroots))
                dblk(i,j,1) = dblk(i,j,1) + der(1)*fac
                dblk(i,j,2) = dblk(i,j,2) + der(2)*fac
                dblk(i,j,3) = dblk(i,j,3) + der(3)*fac
            END DO
        END DO
        END ASSOCIATE
    END DO
 END SUBROUTINE

!> @brief Charge-center (Hellmann-Feynman) first derivative of the
!>   nuclear-attraction integral for one charge, returned as an (inao, jnao, 3)
!>   block (not contracted). Mirrors comp_coulomb_helfeyder1.
 SUBROUTINE comp_coulomb_helfeyder1_block(cp, c, znuc, dblk)
    TYPE(shpair_t), INTENT(IN) :: cp
    REAL(REAL64), INTENT(IN) :: c(3), znuc
    REAL(REAL64), CONTIGUOUS, INTENT(INOUT) :: dblk(:,:,:)

    REAL(REAL64) :: xx, fac, der(3), ric(3)
    type(rys_root_t) :: ryscomp
    INTEGER :: id, i, j, ix, iy, iz, jx, jy, jz
    real(real64) :: xyzin(0:2*max_ang+1, 0:max_ang+1,3,max_nroots)
    real(real64) :: dxyzc(0:max_ang_pad,0:max_ang,3,max_nroots)

    ric = cp%ri(:3) - c(:3)
    DO id = 1, cp%numpairs
        ASSOCIATE (pp => cp%p(id), iang => cp%iang, jang => cp%jang, &
                   inao => cp%inao, jnao => cp%jnao)
        xx = pp%aa*sum((pp%r-c)**2)
        fac = pp%expfac*TWOPI*2
        ryscomp%nroots = cp%nroots
        ryscomp%x = xx
        CALL DQGaussRys(ryscomp, cp, id, c, znuc, xyzin)
        CALL der_helfey_xyz(dxyzc, xyzin, iang, jang, ric, cp%nroots)
        DO i = 1, inao
            ix = CART_X(i,iang); iy = CART_Y(i,iang); iz = CART_Z(i,iang)
            DO j = 1, jnao
                jx = CART_X(j,jang); jy = CART_Y(j,jang); jz = CART_Z(j,jang)
                der(1) = sum(dxyzc(jx,ix,1,1:cp%nroots)*xyzin(jy,iy,2,1:cp%nroots)*xyzin(jz,iz,3,1:cp%nroots))
                der(2) = sum(xyzin(jx,ix,1,1:cp%nroots)*dxyzc(jy,iy,2,1:cp%nroots)*xyzin(jz,iz,3,1:cp%nroots))
                der(3) = sum(xyzin(jx,ix,1,1:cp%nroots)*xyzin(jy,iy,2,1:cp%nroots)*dxyzc(jz,iz,3,1:cp%nroots))
                dblk(i,j,1) = dblk(i,j,1) + der(1)*fac
                dblk(i,j,2) = dblk(i,j,2) + der(2)*fac
                dblk(i,j,3) = dblk(i,j,3) + der(3)*fac
            END DO
        END DO
        END ASSOCIATE
    END DO
 END SUBROUTINE

!> @brief Uncontracted per-AO basis-center second-derivative blocks of the
!>   nuclear-attraction integral for one charge center c (Gate 2 shared kernel).
!> @details For a shell pair (bra X on atom A, ket Y on atom B) and charge centre
!>   c this returns, for every cartesian AO pair (i in bra, j in ket):
!>     pAA(a,b,i,j) = d2/dA_a dA_b <i| znuc/|r-c| |j>   (bra-bra, symmetric in a,b)
!>     pAB(a,b,i,j) = d2/dA_a dB_b <i| znuc/|r-c| |j>   (bra-ket mixed, NOT symmetric)
!>
!>   The production basis-basis second derivative uses angular-momentum (AM) shift
!>   identities and therefore does NOT differentiate the Rys roots/weights: the
!>   bra second derivative is der2_coul_xyz (the bra raise/lower recursion applied
!>   twice), the ket first derivative and the mixed bra-ket second derivative are
!>   the analogous explicit index recurrences. The roots/weights are merely
!>   recomputed for the shifted integral class at the corrected second-derivative
!>   count nroots_der2 = floor((Li+Lj+2)/2)+1 (set here; NOT inherited from
!>   cp%nroots). The validated rys_deriv.F90 layer is not used on this path.
!>
!>   Blocks are returned in the unnormalized cartesian convention (apply the
!>   bfnrm shell normalization at contraction time, exactly as comp_coulomb_der1
!>   does for the gradient). Only s/p/d/f shells (nroots_der2 <= 5, i.e. the
!>   closed-form rys_rt1..rys_rt5 regime) are supported; larger shells abort.
 SUBROUTINE comp_coulomb_der2_blocks(cp, c, znuc, pAA, pAB)
    TYPE(shpair_t), INTENT(IN) :: cp
    REAL(REAL64), INTENT(IN) :: c(3), znuc
    REAL(REAL64), INTENT(OUT) :: pAA(:,:,:,:), pAB(:,:,:,:)   ! (3,3,inao,jnao)

    INTEGER :: id, i, j, nr, ix, iy, iz, jx, jy, jz
    INTEGER :: nroots_der2
    REAL(REAL64) :: xx, fac, aj
    type(rys_root_t) :: ryscomp
    real(real64) :: xyzin(0:2*max_ang+2, 0:max_ang+2, 3, max_nroots)
    real(real64) :: gDA (0:max_ang+2, 0:max_ang+2, 3, max_nroots)
    real(real64) :: gDAA(0:max_ang+2, 0:max_ang+2, 3, max_nroots)
    real(real64) :: dket(3, max_nroots)   ! per-root ket (B-center) 1D first derivatives
    real(real64) :: d2ab(3, max_nroots)   ! per-root mixed bra-ket 1D second derivatives

    nroots_der2 = (cp%iang + cp%jang + 2)/2 + 1
    if (nroots_der2 > 5) &
        error stop 'comp_coulomb_der2_blocks: shell L>=4 not supported (nroots_der2>5)'

    pAA = 0.0_real64
    pAB = 0.0_real64

    DO id = 1, cp%numpairs
        ASSOCIATE (pp => cp%p(id), &
                   iang => cp%iang, jang => cp%jang, &
                   inao => cp%inao, jnao => cp%jnao)

        xx = pp%aa * sum((pp%r - c)**2)
        ryscomp%nroots = nroots_der2
        ryscomp%x = xx
        CALL QGaussRys(ryscomp, cp, id, c, znuc, xyzin, 2)

        ! Bra first and second derivatives at the correct quadrature order
        CALL der_coul_xyz (gDA,  xyzin, iang, jang, pp%ai, nroots_der2)
        CALL der2_coul_xyz(gDAA, xyzin, iang, jang, pp%ai, nroots_der2)

        fac = pp%expfac * TWOPI * pp%aa1
        aj = pp%aj
        nr = nroots_der2

        DO i = 1, inao
            ix = CART_X(i,iang); iy = CART_Y(i,iang); iz = CART_Z(i,iang)
            DO j = 1, jnao
                jx = CART_X(j,jang); jy = CART_Y(j,jang); jz = CART_Z(j,jang)

                ! Ket (B-center) first derivative per root: d/dB_q = 2*aj*[j+1,...] - j*[j-1,...]
                dket(1,1:nr) = 2*aj*xyzin(jx+1,ix,1,1:nr)
                if (jx > 0) dket(1,1:nr) = dket(1,1:nr) - jx*xyzin(jx-1,ix,1,1:nr)
                dket(2,1:nr) = 2*aj*xyzin(jy+1,iy,2,1:nr)
                if (jy > 0) dket(2,1:nr) = dket(2,1:nr) - jy*xyzin(jy-1,iy,2,1:nr)
                dket(3,1:nr) = 2*aj*xyzin(jz+1,iz,3,1:nr)
                if (jz > 0) dket(3,1:nr) = dket(3,1:nr) - jz*xyzin(jz-1,iz,3,1:nr)

                ! Mixed bra-ket second derivative per root:
                ! d2/dA_q dB_q = 2*ai*(2*aj*f(j+1,i+1)-j*f(j-1,i+1)) - i*(2*aj*f(j+1,i-1)-j*f(j-1,i-1))
                d2ab(1,1:nr) = 2*pp%ai * 2*aj * xyzin(jx+1,ix+1,1,1:nr)
                if (jx > 0)          d2ab(1,1:nr) = d2ab(1,1:nr) - 2*pp%ai*jx*xyzin(jx-1,ix+1,1,1:nr)
                if (ix > 0)          d2ab(1,1:nr) = d2ab(1,1:nr) - ix*2*aj*xyzin(jx+1,ix-1,1,1:nr)
                if (ix > 0 .and. jx > 0) d2ab(1,1:nr) = d2ab(1,1:nr) + ix*jx*xyzin(jx-1,ix-1,1,1:nr)

                d2ab(2,1:nr) = 2*pp%ai * 2*aj * xyzin(jy+1,iy+1,2,1:nr)
                if (jy > 0)          d2ab(2,1:nr) = d2ab(2,1:nr) - 2*pp%ai*jy*xyzin(jy-1,iy+1,2,1:nr)
                if (iy > 0)          d2ab(2,1:nr) = d2ab(2,1:nr) - iy*2*aj*xyzin(jy+1,iy-1,2,1:nr)
                if (iy > 0 .and. jy > 0) d2ab(2,1:nr) = d2ab(2,1:nr) + iy*jy*xyzin(jy-1,iy-1,2,1:nr)

                d2ab(3,1:nr) = 2*pp%ai * 2*aj * xyzin(jz+1,iz+1,3,1:nr)
                if (jz > 0)          d2ab(3,1:nr) = d2ab(3,1:nr) - 2*pp%ai*jz*xyzin(jz-1,iz+1,3,1:nr)
                if (iz > 0)          d2ab(3,1:nr) = d2ab(3,1:nr) - iz*2*aj*xyzin(jz+1,iz-1,3,1:nr)
                if (iz > 0 .and. jz > 0) d2ab(3,1:nr) = d2ab(3,1:nr) + iz*jz*xyzin(jz-1,iz-1,3,1:nr)

                ! p_AA = d2/dA^2 (symmetric; fill then mirror a<->b)
                pAA(1,1,i,j) = pAA(1,1,i,j) + fac*sum(gDAA(jx,ix,1,1:nr)*xyzin(jy,iy,2,1:nr)*xyzin(jz,iz,3,1:nr))
                pAA(2,2,i,j) = pAA(2,2,i,j) + fac*sum(xyzin(jx,ix,1,1:nr)*gDAA(jy,iy,2,1:nr)*xyzin(jz,iz,3,1:nr))
                pAA(3,3,i,j) = pAA(3,3,i,j) + fac*sum(xyzin(jx,ix,1,1:nr)*xyzin(jy,iy,2,1:nr)*gDAA(jz,iz,3,1:nr))
                pAA(2,1,i,j) = pAA(2,1,i,j) + fac*sum(gDA(jx,ix,1,1:nr)*gDA(jy,iy,2,1:nr)*xyzin(jz,iz,3,1:nr))
                pAA(3,1,i,j) = pAA(3,1,i,j) + fac*sum(gDA(jx,ix,1,1:nr)*xyzin(jy,iy,2,1:nr)*gDA(jz,iz,3,1:nr))
                pAA(3,2,i,j) = pAA(3,2,i,j) + fac*sum(xyzin(jx,ix,1,1:nr)*gDA(jy,iy,2,1:nr)*gDA(jz,iz,3,1:nr))
                pAA(1,2,i,j) = pAA(2,1,i,j)
                pAA(1,3,i,j) = pAA(3,1,i,j)
                pAA(2,3,i,j) = pAA(3,2,i,j)

                ! p_AB = d2/dA_a dB_b (NOT symmetric; pAB(a,b) = d2/dA_a dB_b)
                pAB(1,1,i,j) = pAB(1,1,i,j) + fac*sum(d2ab(1,1:nr)*xyzin(jy,iy,2,1:nr)*xyzin(jz,iz,3,1:nr))
                pAB(2,2,i,j) = pAB(2,2,i,j) + fac*sum(xyzin(jx,ix,1,1:nr)*d2ab(2,1:nr)*xyzin(jz,iz,3,1:nr))
                pAB(3,3,i,j) = pAB(3,3,i,j) + fac*sum(xyzin(jx,ix,1,1:nr)*xyzin(jy,iy,2,1:nr)*d2ab(3,1:nr))
                pAB(1,2,i,j) = pAB(1,2,i,j) + fac*sum(gDA(jx,ix,1,1:nr)*dket(2,1:nr)*xyzin(jz,iz,3,1:nr))
                pAB(2,1,i,j) = pAB(2,1,i,j) + fac*sum(dket(1,1:nr)*gDA(jy,iy,2,1:nr)*xyzin(jz,iz,3,1:nr))
                pAB(1,3,i,j) = pAB(1,3,i,j) + fac*sum(gDA(jx,ix,1,1:nr)*xyzin(jy,iy,2,1:nr)*dket(3,1:nr))
                pAB(3,1,i,j) = pAB(3,1,i,j) + fac*sum(dket(1,1:nr)*xyzin(jy,iy,2,1:nr)*gDA(jz,iz,3,1:nr))
                pAB(2,3,i,j) = pAB(2,3,i,j) + fac*sum(xyzin(jx,ix,1,1:nr)*gDA(jy,iy,2,1:nr)*dket(3,1:nr))
                pAB(3,2,i,j) = pAB(3,2,i,j) + fac*sum(xyzin(jx,ix,1,1:nr)*dket(2,1:nr)*gDA(jz,iz,3,1:nr))

            END DO
        END DO
        END ASSOCIATE
    END DO

  END SUBROUTINE

!> @brief Bra-bra and bra-charge second derivatives of the nuclear-attraction
!>   integral for a single charge center c, contracted with a density block.
!> @details Thin contraction wrapper over comp_coulomb_der2_blocks (the shared
!>   production AM-shift kernel). On output:
!>     p_XX += sum_ij dij(i,j) * d2/dX^2 <i|znuc/|r-c||j>        (bra-bra)
!>     p_XC += -sum_ij dij(i,j) * (d2/dX^2 + d2/dX dY) <i|...|j> (bra-charge)
!>   The bra-charge block uses single-center translational invariance
!>   d/dc = -(d/dX + d/dY); no Rys root/weight differentiation is involved
!>   (see comp_coulomb_der2_blocks). The caller (hess_en) obtains p_YY, p_YC by
!>   calling again with the shell pair swapped, then assembles the 9 atom blocks.
 SUBROUTINE comp_coulomb_der2_braC(cp, c, znuc, dij, p_XX, p_XC)
    TYPE(shpair_t), INTENT(IN) :: cp
    REAL(REAL64), INTENT(IN) :: c(3), znuc
    REAL(REAL64), INTENT(IN) :: dij(:,:)
    REAL(REAL64), CONTIGUOUS, INTENT(INOUT) :: p_XX(:,:), p_XC(:,:)

    REAL(REAL64) :: pAA(3,3,cp%inao,cp%jnao), pAB(3,3,cp%inao,cp%jnao)
    REAL(REAL64) :: bAA(3,3), bAB(3,3)
    INTEGER :: i, j

    CALL comp_coulomb_der2_blocks(cp, c, znuc, pAA, pAB)

    bAA = 0.0_real64
    bAB = 0.0_real64
    DO i = 1, cp%inao
        DO j = 1, cp%jnao
            bAA = bAA + dij(i,j)*pAA(:,:,i,j)
            bAB = bAB + dij(i,j)*pAB(:,:,i,j)
        END DO
    END DO

    p_XX = p_XX + bAA
    ! p_XC = d2/dA dC = -(d2/dA^2 + d2/dA dB) by translational invariance
    p_XC = p_XC - (bAA + bAB)

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

!> @brief Second derivative of the 1D Coulomb (nuclear-attraction) integrals
!>  with respect to the bra center, obtained by applying the bra-center
!>  derivative recursion (der_coul_xyz) twice:
!>    d2[j,i] = 4 ai^2 [j,i+2] - 2 ai (2i+1) [j,i] + i(i-1) [j,i-2]
!>  per Rys root. The input array must be available up to bra index lit+2
!>  (build QGaussRys with igrd=2 and a correspondingly sized xyzin).
!>
!>  Together with the analogous ket-center derivatives, this provides the
!>  basis-center second-derivative blocks (AA, AB, BB) of the nuclear-attraction
!>  Hessian. The charge-center (Hellmann-Feynman) and mixed blocks follow from
!>  translational invariance, d/dC = -(d/dA + d/dB), so no second-derivative Rys
!>  root machinery is required.
 SUBROUTINE der2_coul_xyz(d2xyz,xyzin,lit,ljt,ai,nroots)
!dir$ attributes forceinline :: der2_coul_xyz
    REAL(REAL64), INTENT(IN) ::  ai
    REAL(REAL64), CONTIGUOUS, INTENT(IN) ::  xyzin(0:,0:,:,:)
    REAL(REAL64), CONTIGUOUS, INTENT(OUT) :: d2xyz(0:,0:,:,:)
    INTEGER, INTENT(IN) :: lit, ljt, nroots

    INTEGER :: i
!dir$ assume_aligned xyzin : 64
!dir$ assume_aligned d2xyz : 64

    d2xyz(0:ljt,0:lit,1:3,1:nroots) = 4*ai*ai * xyzin(0:ljt,2:lit+2,1:3,1:nroots)

    DO i = 0, lit
        d2xyz(0:ljt,i,1:3,1:nroots) = d2xyz(0:ljt,i,1:3,1:nroots) &
            - 2*ai*(2*i+1)*xyzin(0:ljt,i,1:3,1:nroots)
    END DO

    DO i = 2, lit
        d2xyz(0:ljt,i,1:3,1:nroots) = d2xyz(0:ljt,i,1:3,1:nroots) &
            + i*(i-1)*xyzin(0:ljt,i-2,1:3,1:nroots)
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

!> @brief Second derivative of 1D overlap/kinetic integrals w.r.t. the bra
!>  center, obtained by applying the bra-center derivative operator twice:
!>    d2[j,i] = 4 ai^2 [j,i+2] - 2 ai (2i+1) [j,i] + i(i-1) [j,i-2]
!>  This is identical to der_kinovl_xyz composed with itself; the input array
!>  must therefore be available up to bra index lit+2.
 SUBROUTINE der2_kinovl_xyz(d2xyz,xyz,lit,ljt,ai)
!dir$ attributes forceinline :: der2_kinovl_xyz
    REAL(REAL64), INTENT(IN) ::  ai
    REAL(REAL64), CONTIGUOUS, INTENT(IN) ::  xyz(0:,0:,:)
    REAL(REAL64), CONTIGUOUS, INTENT(OUT) :: d2xyz(0:,0:,:)
    INTEGER, INTENT(IN) :: lit, ljt
    INTEGER :: i

    d2xyz(0:ljt,0:lit,:) = 4*ai*ai * xyz(0:ljt,2:lit+2,:)

    DO i = 0, lit
        d2xyz(0:ljt,i,:) = d2xyz(0:ljt,i,:) - 2*ai*(2*i+1)*xyz(0:ljt,i,:)
    END DO

    DO i = 2, lit
        d2xyz(0:ljt,i,:) = d2xyz(0:ljt,i,:) + i*(i-1)*xyz(0:ljt,i-2,:)
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
 SUBROUTINE DQGaussRys(ryscomp, cp, id, c, znuc, xyzin, igrd)
!dir$ attributes forceinline :: DQGaussRys
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

    igrd1 = 1
    IF (present(igrd)) igrd1 = igrd

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
        DO nj = 2, iang+jang+igrd1
            xyzin(nj,0,:,k) = d(:)*xyzin(nj-1,0,:,k) + (nj-1)*b*xyzin(nj-2,0,:,k)
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

!> @brief Compute double derivative of 1D Coulomb integral table
!>        for PVP integrals: d²/d(bra_x) d(ket_x) of xyzin
!>
!> @details For each direction gamma and each Rys root k:
!>
!>   dxyz(j,i,gamma,k) =
!>       4*ai*aj * xyzin(j+1, i+1, gamma, k)
!>     - 2*ai*j  * xyzin(j-1, i+1, gamma, k)   ! j=0 => 0
!>     - 2*aj*i  * xyzin(j+1, i-1, gamma, k)   ! i=0 => 0
!>     +     i*j * xyzin(j-1, i-1, gamma, k)   ! i=0 or j=0 => 0
!>
!>   where i = 0..li  (bra angular momentum)
!>         j = 0..lj  (ket angular momentum)
!>
!> @param[in]  xyzin   1D Coulomb integral table, built with igrd=2
!> @param[in]  li      bra max angular momentum
!> @param[in]  lj      ket max angular momentum
!> @param[in]  ai      bra primitive exponent
!> @param[in]  aj      ket primitive exponent
!> @param[in]  nroots  number of Rys roots
!> @param[out] dxyz    derivative table, same indexing as xyzin
!
!> @author   Vladimir Makhnev
!
 subroutine pvp_xyz_ij(xyzin, li, lj, ai, aj, nroots, dxyz)
!dir$ attributes forceinline :: pvp_xyz_ij

    implicit none

    real(real64), contiguous, intent(in)  :: xyzin(0:,0:,:,:)
    real(real64), contiguous, intent(out) :: dxyz(0:,0:,:,:)
    real(real64), intent(in) :: ai, aj
    integer,      intent(in) :: li, lj, nroots

    integer      :: i, j, k
    real(real64) :: ai2, aj2

    ai2 = 2.0_real64 * ai
    aj2 = 2.0_real64 * aj

    do k = 1, nroots

        ! i=0, j=0: only 4*ai*aj term survives
        dxyz(0, 0, :, k) = ai2*aj2 * xyzin(1, 1, :, k)

        ! i=0, j>0: j-1 terms vanish 
        ! only 4ai*aj and -2ai*j terms survive
        do j = 1, lj
            dxyz(j, 0, :, k) =  ai2*aj2 * xyzin(j+1, 1, :, k) &
                               - ai2*j   * xyzin(j-1, 1, :, k)
        end do

        ! j=0, i>0: i-1 terms vanish
        ! only 4ai*aj and -2aj*i terms survive
        do i = 1, li
            dxyz(0, i, :, k) =  ai2*aj2 * xyzin(1, i+1, :, k) &
                               - aj2*i   * xyzin(1, i-1, :, k)
        end do

        ! general case i>0, j>0: all four terms
        do i = 1, li
            do j = 1, lj
                dxyz(j, i, :, k) =  ai2*aj2 * xyzin(j+1, i+1, :, k) &
                                   - ai2*j   * xyzin(j-1, i+1, :, k) &
                                   - aj2*i   * xyzin(j+1, i-1, :, k) &
                                   +  real(i*j, real64) * xyzin(j-1, i-1, :, k)
            end do
        end do

    end do

 end subroutine pvp_xyz_ij


!> @brief Compute primitive block of PVP integrals
!>        <mu | p . (-Z/|r-C|) . p | nu>
!>        = <d(mu)/dx | -Z/|r-C| | d(nu)/dx>
!>        + <d(mu)/dy | -Z/|r-C| | d(nu)/dy>
!>        + <d(mu)/dz | -Z/|r-C| | d(nu)/dz>
!>
!> @param[in]     cp     shell pair data
!> @param[in]     id     current pair of primitives
!> @param[in]     c      coordinates of the nucleus
!> @param[in]     znuc   nuclear charge (passed as -Z, same as coulomb)
!> @param[inout]  pvpblk block of PVP integrals (accumulated)
!
!> @author  Vladimir Makhnev
!
 subroutine comp_pvp_int1_prim(cp, id, c, znuc, pvpblk)
!dir$ attributes inline :: comp_pvp_int1_prim

    implicit none

    type(shpair_t),  intent(in)    :: cp
    integer,         intent(in)    :: id
    real(real64),    intent(in)    :: c(3), znuc
    real(real64), contiguous, intent(inout) :: pvpblk(:)
!dir$ assume_aligned pvpblk : 64
    real(real64) :: xyzin(0:2*max_ang+3, 0:max_ang+2, 3, max_nroots+1)
    real(real64) :: dxyz (0:2*max_ang+3, 0:max_ang+2, 3, max_nroots+1)
!dir$ assume_aligned xyzin : 64
!dir$ assume_aligned dxyz  : 64

    type(rys_root_t) :: ryscomp
    integer      :: i, j, ij, jmax
    integer      :: nx, ny, nz, mx, my, mz
    integer      :: nroots_pvp
    real(real64) :: xx, dij, dum

    associate( pp   => cp%p(id),         &
               iang => cp%iang, jang => cp%jang, &
               inao => cp%inao, jnao => cp%jnao  )

    nroots_pvp = cp%nroots + 1

    xx = pp%aa * sum((pp%r - c)**2)

    ryscomp%nroots = nroots_pvp
    ryscomp%x      = xx
    call QGaussRys(ryscomp, cp, id, c, znuc, xyzin, 2)

    call pvp_xyz_ij(xyzin, iang, jang, pp%ai, pp%aj, nroots_pvp, dxyz)

    dij = pp%expfac * TWOPI * pp%aa1

    ij   = 0
    jmax = jnao
    do i = 1, inao
        nx = CART_X(i, iang)
        ny = CART_Y(i, iang)
        nz = CART_Z(i, iang)
        if (cp%iandj) jmax = i
        do j = 1, jmax
            mx = CART_X(j, jang)
            my = CART_Y(j, jang)
            mz = CART_Z(j, jang)
            ij = ij + 1

            ! (pVp) = p_x V p_x  +  p_y V p_y  +  p_z V p_z
            dum = sum(  dxyz (mx, nx, 1, 1:nroots_pvp)   &  ! d²/dx_bra dx_ket
                      * xyzin(my, ny, 2, 1:nroots_pvp)   &  ! y: 
                      * xyzin(mz, nz, 3, 1:nroots_pvp) ) &  ! z:
                + sum(  xyzin(mx, nx, 1, 1:nroots_pvp)   &  ! x: 
                      * dxyz (my, ny, 2, 1:nroots_pvp)   &  ! d²/dy_bra dy_ket
                      * xyzin(mz, nz, 3, 1:nroots_pvp) ) &  ! z:
                + sum(  xyzin(mx, nx, 1, 1:nroots_pvp)   &  ! x:
                      * xyzin(my, ny, 2, 1:nroots_pvp)   &  ! y:
                      * dxyz (mz, nz, 3, 1:nroots_pvp) )    ! d²/dz_bra dz_ket

            pvpblk(ij) = pvpblk(ij) + dij * dum

        end do
    end do

    end associate

 end subroutine comp_pvp_int1_prim

!> @brief Build one-sided Gaussian derivative tables for SOC integrals
!> @details
!>   di(m,n) = 2*ai*xyzin(m,n+1) - n*xyzin(m,n-1)  [d/d bra, second index]
!>   dj(m,n) = 2*aj*xyzin(m+1,n) - m*xyzin(m-1,n)  [d/d ket, first index]
!> @param[in]  xyzin   1D Rys integral table (dims: (Lj,Li,XYZ,NRoots))
!> @param[in]  iang    angular momentum of bra shell
!> @param[in]  jang    angular momentum of ket shell
!> @param[in]  ai      exponent of bra primitive
!> @param[in]  aj      exponent of ket primitive
!> @param[in]  nroots  number of Rys roots
!> @param[out] di      bra-side derivative table (dims: (Lj,Li,XYZ,NRoots))
!> @param[out] dj      ket-side derivative table (dims: (Lj,Li,XYZ,NRoots))
!
!> @author  Vladimir Makhnev
!> @date    March 2026

 subroutine soc_xyz_ij(xyzin, iang, jang, ai, aj, nroots, di, dj)
!dir$ attributes forceinline :: soc_xyz_ij
    implicit none
    real(real64), contiguous, intent(in)  :: xyzin(0:, 0:, :, :)
    real(real64), contiguous, intent(out) :: di(0:, 0:, :, :)
    real(real64), contiguous, intent(out) :: dj(0:, 0:, :, :)
    integer,      intent(in) :: iang, jang, nroots
    real(real64), intent(in) :: ai, aj

    integer :: n
!dir$ assume_aligned xyzin : 64
!dir$ assume_aligned di    : 64
!dir$ assume_aligned dj    : 64

! bra derivative (second index): di(m,n) = n*xyzin(m,n-1) - 2*ai*xyzin(m,n+1)
    di(0:jang, 0:iang, 1:3, 1:nroots) = -2*ai * xyzin(0:jang, 1:iang+1, 1:3, 1:nroots)
    do n = 1, iang
        di(0:jang, n, 1:3, 1:nroots) = di(0:jang, n, 1:3, 1:nroots) &
                                      + n * xyzin(0:jang, n-1, 1:3, 1:nroots)
    end do

    ! ket derivative (first index): dj(m,n) = m*xyzin(m-1,n) - 2*aj*xyzin(m+1,n)
    dj(0:jang, 0:iang, 1:3, 1:nroots) = -2*aj * xyzin(1:jang+1, 0:iang, 1:3, 1:nroots)
    do n = 1, jang
        dj(n, 0:iang, 1:3, 1:nroots) = dj(n, 0:iang, 1:3, 1:nroots) &
                                      + n * xyzin(n-1, 0:iang, 1:3, 1:nroots)
    end do

 end subroutine soc_xyz_ij


!> @brief Compute primitive block of 1e SOC integrals for one nucleus
!> @details Evaluates <mu|Z_eff*L/r_A^3|nu> via IBP reduction to Coulomb
!>  integrals with one-sided Gaussian derivatives:
!>   Lx: <d/dy mu|1/r|d/dz nu> - <d/dz mu|1/r|d/dy nu>
!>   Ly: <d/dz mu|1/r|d/dx nu> - <d/dx mu|1/r|d/dz nu>
!>   Lz: <d/dx mu|1/r|d/dy nu> - <d/dy mu|1/r|d/dx nu>
!>
!> @param[in]     cp      shell pair data
!> @param[in]     id      current pair of primitives
!> @param[in]     c       nuclear coordinates
!> @param[in]     znuc    effective nuclear charge Z_eff
!> @param[inout]  socblk  accumulated SOC block (nfunc,3): Lx, Ly, Lz
!
!> @author  Vladimir Makhnev
!> @date    March 2026

 subroutine comp_soc_int1_prim(cp, id, c, znuc, socblk)
!dir$ attributes inline :: comp_soc_int1_prim
    implicit none

    type(shpair_t), intent(in)              :: cp
    integer,        intent(in)              :: id
    real(real64),   intent(in)              :: c(3), znuc
    real(real64), contiguous, intent(inout) :: socblk(:, :)

    type(rys_root_t) :: ryscomp
    integer      :: i, j, ij, jmax, nroots_soc
    integer      :: nx, ny, nz, mx, my, mz
    real(real64) :: dij, lx, ly, lz

    real(real64) :: xyzin(0:2*max_ang+1, 0:max_ang+1, 3, max_nroots+1)
    real(real64) :: di(0:max_ang_pad, 0:max_ang, 3, max_nroots+1)
    real(real64) :: dj(0:max_ang_pad, 0:max_ang, 3, max_nroots+1)
!dir$ assume_aligned xyzin  : 64
!dir$ assume_aligned socblk : 64

    associate (pp   => cp%p(id), &
               iang => cp%iang,  jang => cp%jang, &
               inao => cp%inao,  jnao => cp%jnao)
    ! 
    nroots_soc     = cp%nroots + 2   
    ryscomp%nroots = nroots_soc
    ryscomp%x      = pp%aa * sum((pp%r - c)**2)
    call QGaussRys(ryscomp, cp, id, c, znuc, xyzin, 2)  


    call soc_xyz_ij(xyzin, iang, jang, pp%ai, pp%aj, nroots_soc, di, dj)

    dij  = pp%expfac * TWOPI * pp%aa1
    ij   = 0
    jmax = jnao

    do i = 1, inao
        nx = CART_X(i, iang);  ny = CART_Y(i, iang);  nz = CART_Z(i, iang)
        if (cp%iandj) jmax = i
        do j = 1, jmax
            mx = CART_X(j, jang);  my = CART_Y(j, jang);  mz = CART_Z(j, jang)
            ij = ij + 1

            ! Lx = <d/dy mu|1/r|d/dz nu> - <d/dz mu|1/r|d/dy nu>
            lx = sum( xyzin(mx,nx,1,1:nroots_soc) &
                    *    di(my,ny,2,1:nroots_soc)  &
                    *    dj(mz,nz,3,1:nroots_soc) )&
               - sum( xyzin(mx,nx,1,1:nroots_soc) &
                    *    di(mz,nz,3,1:nroots_soc)  &
                    *    dj(my,ny,2,1:nroots_soc) )

            ! Ly = <d/dz mu|1/r|d/dx nu> - <d/dx mu|1/r|d/dz nu>
            ly = sum( xyzin(my,ny,2,1:nroots_soc) &
                    *    di(mz,nz,3,1:nroots_soc)  &
                    *    dj(mx,nx,1,1:nroots_soc) )&
               - sum( xyzin(my,ny,2,1:nroots_soc) &
                    *    di(mx,nx,1,1:nroots_soc)  &
                    *    dj(mz,nz,3,1:nroots_soc) )

            ! Lz = <d/dx mu|1/r|d/dy nu> - <d/dy mu|1/r|d/dx nu>
            lz = sum( xyzin(mz,nz,3,1:nroots_soc) &
                    *    di(mx,nx,1,1:nroots_soc)  &
                    *    dj(my,ny,2,1:nroots_soc) )&
               - sum( xyzin(mz,nz,3,1:nroots_soc) &
                    *    di(my,ny,2,1:nroots_soc)  &
                    *    dj(mx,nx,1,1:nroots_soc) )

            socblk(ij, 1) = socblk(ij, 1) + dij * lx
            socblk(ij, 2) = socblk(ij, 2) + dij * ly
            socblk(ij, 3) = socblk(ij, 3) + dij * lz

        end do
    end do

    end associate

 end subroutine comp_soc_int1_prim

! 2e part starts here. 

!> @brief Compute 1D integrals for two-electron spin-orbit operator (4-centre ERI)
!> @details Analogue of QGaussRys but for a Gaussian charge distribution
!>   (shell pair cpkl) instead of a point nucleus.  The output gfull is the
!>   4-index table
!>
!>     gfull(nj, ni, nl, nk, xyz, t)
!>       = G_{ni,nj | nk,nl}^{xyz} for Rys root t
!>
!>   after both electron-1 and electron-2 HRR transfers, exactly as
!>   GAMESS XYZ2E builds XINT(1+NI+MAXP1*NJ, 1+NK+MAXP*NL).
!>
!>   Calling convention in comp_soc_int2_prim:
!>     for each (k,l) subshell with angular powers (nxk,nyk,nzk), (nxl,nyl,nzl):
!>       xyzin(nj,ni,1,t) = gfull(nj,ni,nxl,nxk,1,t)
!>       xyzin(nj,ni,2,t) = gfull(nj,ni,nyl,nyk,2,t)
!>       xyzin(nj,ni,3,t) = gfull(nj,ni,nzl,nzk,3,t)
!>     call soc_xyz_ij(xyzin, ...) as in the 1e case.
!>
!>   Corresponds to GAMESS SOINT2 + XYZ2E.
!>
!> @param[inout] ryscomp   Rys object; caller sets nroots; x is set here
!> @param[in]    cpij      shell pair electron 1 (bra=I, ket=J)
!> @param[in]    idij      primitive index in cpij
!> @param[in]    cpkl      shell pair electron 2 (bra=K, ket=L)
!> @param[in]    idkl      primitive index in cpkl
!> @param[out]   gfull     4-index integral table
!>                         shape (0:jang+1, 0:iang+1, 0:lang, 0:kang, 3, nroots)

 SUBROUTINE QGaussRys2e(ryscomp, cpij, idij, cpkl, idkl, gfull)
!dir$ attributes forceinline :: QGaussRys2e
    USE ISO_FORTRAN_ENV, ONLY: real64
    USE mod_shell_tools, ONLY: shpair_t
    USE rys,             ONLY: rys_root_t
    IMPLICIT NONE

    TYPE(rys_root_t), INTENT(INOUT) :: ryscomp
    TYPE(shpair_t),   INTENT(IN)    :: cpij, cpkl
    INTEGER,          INTENT(IN)    :: idij, idkl
    REAL(real64), CONTIGUOUS, INTENT(OUT) :: gfull(0:, 0:, 0:, 0:, :, :)
!dir$ assume_aligned gfull : 64
    ! --- local constants ---
    REAL(real64), PARAMETER :: PI252 = 34.986836655250_real64   ! 2*pi^(5/2)

    ! --- local scalars ---
    INTEGER      :: t, n, m, ni, nj, nk, nl, nmax
    REAL(real64) :: f00, aandb, rho
    REAL(real64) :: b00, b10, bp01, c10, cp01, cp10, c01
    REAL(real64) :: expe
    REAL(real64) :: c00(3), cp00(3), dij(3), dkl(3), pq(3)
    integer :: nj_e1, ni_max, n_lo, n_hi
    ! --- intermediate 2D table (GAMESS-style packed: row=electron1, col=electron2) ---
    ! Dimensions: (0:maxij+1, 0:maxkl+1, 3)
    ! Row index N = NI + maxp1*NJ  where maxp1 = maxij+1
    ! Col index M = NK + maxp *NL  where maxp  = maxkl+1
    INTEGER      :: maxij, maxkl, maxp1, maxp
    INTEGER, PARAMETER :: ND52 = (MAX_ANG+2)**2   ! electron 1, row
    INTEGER, PARAMETER :: ND51 = (MAX_ANG+1)**2   ! electron 2, col
    REAL(real64) :: g(0:ND52, 0:ND51, 3)
    ASSOCIATE ( ppij => cpij%p(idij), ppkl => cpkl%p(idkl), &
                iang => cpij%iang,    jang => cpij%jang,     &
                kang => cpkl%iang,    lang => cpkl%jang )

    ! ----------------------------------------------------------------
    ! Geometry
    ! ----------------------------------------------------------------
    aandb = ppij%aa + ppkl%aa
    rho   = ppij%aa * ppkl%aa / aandb
    pq    = ppij%r  - ppkl%r
    ryscomp%x = rho * SUM(pq**2)

    CALL ryscomp%evaluate()

    maxij = iang + jang + 1   ! = MAXIJ in GAMESS
    maxkl = kang + lang        ! = MAXKL in GAMESS
    maxp1 = maxij + 1 ! MAX_ANG + 2 !maxij + 1          ! stride for NJ in row index
    maxp  = maxkl + 1 ! MAX_ANG + 1 !maxkl + 1          ! stride for NL in col index


!    print *, 'maxij=', maxij, ' maxp1=', maxp1, ' g size dim1=', maxij+2
    ! Shifts for HRR (GAMESS: DXIJ, DXKL etc.)
    dij = cpij%ri - cpij%rj   ! A - B
    dkl = cpkl%ri - cpkl%rj   ! C - D

    ! Prefactor without Rys weight (GAMESS: EXPE without W(t))
    expe = PI252 / (ppij%aa * ppkl%aa * SQRT(aandb)) &
         * EXP(-ryscomp%x / rho)!* ppij%expfac * ppkl%expfac
!write(*,'(a,5e20.12)') 'OQP expe rho x ai aj:', &
!    expe, rho, ryscomp%x, ppij%ai, ppij%aj
!write(*,'(a,4e20.12)') 'OQP ak al aa bb:', &
!    ppkl%ai, ppkl%aj, ppij%aa, ppkl%aa
!write(*,'(a,2e20.12)') 'OQP Kij Kkl:', ppij%expfac, ppkl%expfac
!    print *, 'DEBUG expe=', expe
    ! ----------------------------------------------------------------
    ! Loop over Rys roots
    ! ----------------------------------------------------------------
    DO t = 1, ryscomp%nroots
        g = 0.0_real64        
        ! F00 = EXPE * W(t)  (GAMESS notation)
        f00 = expe * ryscomp%w(t)
!if (t == 1) write(*,'(a,3e20.12)') 'OQP f00 w(1) t=1:', f00, ryscomp%w(t), ryscomp%u(t)
        ! VRR recurrence coefficients (GAMESS: B00, B10, BP01, XC00, XCP00)
        ! denom = 2*(aa*bb + u*rho*(aa+bb))
        ASSOCIATE (uu => ryscomp%u(t))
        b00  = uu*rho          / (2*(ppij%aa*ppkl%aa + uu*rho*aandb))
        b10  = (ppkl%aa + uu*rho) / (2*(ppij%aa*ppkl%aa + uu*rho*aandb))
        bp01 = (ppij%aa + uu*rho) / (2*(ppij%aa*ppkl%aa + uu*rho*aandb))
        END ASSOCIATE
!        if (t==1) print *, 'b00=', b00, ' b10=', b10, ' bp01=', bp01
        ! VRR centres
        c00  = (ppij%r - cpij%ri) + 2*b00*ppkl%aa * pq   ! XC00
        cp00 = (ppkl%r - cpkl%ri) - 2*b00*ppij%aa * pq   ! XCP00

        ! ----------------------------------------------------------------
        ! Seed values  (GAMESS XYZ2E: XINT(1,1) etc.)
        ! In OQP packed: N=0 → NI=0,NJ=0; M=0 → NK=0,NL=0
        ! z-component carries F00 (the Rys weight × prefactor)
        ! ----------------------------------------------------------------
        g(0, 0, 1) = 1.0_real64
        g(0, 0, 2) = 1.0_real64
        g(0, 0, 3) = f00

        g(1, 0, 1) = c00(1)
        g(1, 0, 2) = c00(2)
        g(1, 0, 3) = c00(3) * f00

        g(0, 1, 1) = cp00(1)
        g(0, 1, 2) = cp00(2)
        g(0, 1, 3) = cp00(3) * f00

        g(1, 1, 1) = c00(1)*cp00(1) + b00
        g(1, 1, 2) = c00(2)*cp00(2) + b00
        g(1, 1, 3) = (c00(3)*cp00(3) + b00) * f00

!        if (t==1) print *, 'DEBUG g(0,0,3)=', g(0,0,3)
!        if (t==1) print *, 'after seed g(1,1,1)=', g(1,1,1)
!        if (t==1) print *, 'after seed g(1,1,3)=', g(1,1,3)
        ! ----------------------------------------------------------------
        ! VRR for electron 1 (N increases, M=0 and M=1)
        ! GAMESS loop 30: C10 = 0; CP10 = B00
        !   G(N+1,0) = C10*G(N-1,0) + C00*G(N,0)   [C10 = (N-1)*B10]
        !   G(N+1,1) = CP10*G(N,0)  + CP00*G(N+1,0) [CP10 = N*B00]
        ! ----------------------------------------------------------------
        c10  = 0.0_real64
        cp10 = b00
        DO n = 2, maxij
            c10  = c10  + b10
            cp10 = cp10 + b00
            g(n, 0, :) = c10*g(n-2, 0, :) + c00*g(n-1, 0, :)
            g(n, 1, :) = cp10*g(n-1, 0, :) + cp00*g(n, 0, :)
        END DO
!        if (t==1) print *, 'after VRR1 g(0,0,3)=', g(0,0,3)
!        if (t==1) print *, 'after VRR1 g(1,1,1)=', g(1,1,1)
        ! ----------------------------------------------------------------
        ! VRR for electron 2 (M increases, N=0 and N=1)
        ! GAMESS loop 60: CP01 = 0; C01 = B00
        !   G(0,M+1) = CP01*G(0,M-1) + CP00*G(0,M)   [CP01 = (M-1)*BP01]
        !   G(1,M+1) = C01*G(0,M)    + C00*G(0,M+1)  [C01  = M*B00]
        ! ----------------------------------------------------------------
        cp01 = 0.0_real64
        c01  = b00
        DO m = 2, maxkl
            cp01 = cp01 + bp01
            c01  = c01  + b00
            g(0, m, :) = cp01*g(0, m-2, :) + cp00*g(0, m-1, :)
            g(1, m, :) = c01*g(0, m-1, :) + c00*g(0, m, :)
            ! mixed recurrence for N >= 2  (GAMESS loop 50)
            ! G(N,M+1) = CP01*G(N,M-1) + CP10_m*G(N-1,M) + CP00*G(N,M)
            ! CP10_m starts at B00 and increments by B00 per N step
            cp10 = b00
            nmax = MIN(maxij, maxij + 2 - m)
            DO n = 2, nmax
                cp10 = cp10 + b00
                g(n, m, :) = cp01*g(n, m-2, :) + cp10*g(n-1, m-1, :) + cp00*g(n, m-1, :)
            END DO
        END DO
!        if (t==1) print *, 'after VRR2 g(0,0,3)=', g(0,0,3)
 
!        if (t==1) print *, 'after VRR2 g(1,1,1)=', g(1,1,1)
        ! ----------------------------------------------------------------
        ! HRR for electron 1: transfer from NI to (NI, NJ) representation
        ! G(NI, NJ) = G(NI+1, NJ-1) + dij * G(NI, NJ-1)
        ! In packed form: g(NI + maxp1*NJ, M) = g(NI + maxp1*(NJ-1) + 1, M)
        !                                       + dij * g(NI + maxp1*(NJ-1), M)
        ! GAMESS: backward loop over NI is required.
        ! ----------------------------------------------------------------
        DO nj = 1, jang + 1
!            if (t==1 .and. nj==1) print *, 'HRR1 start nj=1 g(0,0,3)=', g(0,0,3)
            DO ni = maxij - nj, 0, -1
!                if (t==1) print *, 'HRR1 nj ni=', nj, ni, ' writes to', ni+maxp1*nj
                g(ni + maxp1*nj, 0:maxkl, :) = &
                    g(ni + maxp1*(nj-1) + 1, 0:maxkl, :) &
                  + spread(dij, 1, maxkl+1) * g(ni + maxp1*(nj-1), 0:maxkl, :)
            END DO
        END DO
!        if (t==1) print *, 'after HRR1 g(0,0,3)=', g(0,0,3)
        
 
!        if (t==1) print *, 'after HRR1 g(1,1,1)=', g(1,1,1)
        ! ----------------------------------------------------------------
        ! HRR for electron 2: transfer from NK to (NK, NL) representation
        ! G(NI_packed, NK, NL) = G(NI_packed, NK+1, NL-1)
        !                       + dkl * G(NI_packed, NK, NL-1)
        ! In col packed: g(N, NK + maxp*NL, :) = g(N, NK + maxp*(NL-1)+1, :)
        !                                        + dkl * g(N, NK+maxp*(NL-1), :)
        ! Loop over valid N values (all electron-1 packed indices).
        ! ----------------------------------------------------------------
!        DO nl = 1, lang
!            DO nk = maxkl - nl, 0, -1
!                g(0:maxij, nk + maxp*nl, :) = &
!                    g(0:maxij, nk + maxp*(nl-1) + 1, :) &
!                  + SPREAD(dkl, 1, maxij+1) * g(0:maxij, nk + maxp*(nl-1), :)
!            END DO
!        END DO

! СТАЛО (все NJ строки как в GAMESS):
         DO nl = 1, lang
             DO nk = maxkl - nl, 0, -1
                 DO nj_e1 = 0, jang + 1
                     ni_max = MIN(iang + 1, maxij - nj_e1)
                     IF (ni_max < 0) EXIT
                     n_lo = maxp1 * nj_e1
                     n_hi = ni_max + n_lo
                     g(n_lo:n_hi, nk + maxp*nl, :) = &
                         g(n_lo:n_hi, nk + maxp*(nl-1) + 1, :) &
                       + SPREAD(dkl, 1, ni_max + 1) * g(n_lo:n_hi, nk + maxp*(nl-1), :)
                 END DO
             END DO
         END DO


!        if (t==1) print *, 'after HRR2 g(0,0,3)=', g(0,0,3)

!        if (t==1) print *, 'after HRR2 g(1,1,1)=', g(1,1,1)
        ! ----------------------------------------------------------------
        ! Extract into gfull(nj, ni, nl, nk, xyz, t)
        ! From packed: g(ni + maxp1*nj, nk + maxp*nl, xyz)
        ! ----------------------------------------------------------------
        DO nl = 0, lang
            DO nk = 0, kang
                DO nj = 0, jang + 1
                    DO ni = 0, iang + 1
                        gfull(nj, ni, nl, nk, :, t) = g(ni + maxp1*nj, nk + maxp*nl, :) !g(nj + maxp1*ni, nl + maxp*nk, :) !g(ni + maxp1*nj, nk + maxp*nl, :)
                    END DO
                END DO
            END DO
        END DO
    END DO   ! Rys roots
!write(*,'(a,e20.12)') 'OQP gfull(0,0,0,0,3,1) [F00]:', gfull(0,0,0,0,3,1)
!write(*,'(a,e20.12)') 'OQP gfull(1,0,0,0,1,1) [C00x]:', gfull(1,0,0,0,1,1)
!write(*,'(a,e20.12)') 'OQP gfull(0,0,1,0,1,1) [CP00x]:', gfull(0,0,1,0,1,1)
!rite(*,'(a,e20.12)') 'OQP gfull(0,1,0,1,1,1) [B00x]:', gfull(0,1,0,1,1,1)
    END ASSOCIATE

 END SUBROUTINE QGaussRys2e




!> @brief Two-electron SOC primitive integral over one primitive pair (idij, idkl)
!>
!> Computes the two-electron mean-field spin-orbit coupling contribution
!> for a single pair of contracted primitives:
!>
!>   socblk(ij) += dij_factor * Lx(i,j,k,l)   summed over all (k,l) of electron 2
!>
!> where Lx = <d/dy mu | 1/r12 | d/dz nu> - <d/dz mu | 1/r12 | d/dy nu>
!>
!> Correspondence with GAMESS SOINT2:
!>   - QGaussRys2e  builds gfull  ←→  XYZ2E builds XINT/YINT/ZINT + XINTI/YINTJ etc.
!>   - gfull(nj,ni,nl,nk,xyz,t)  ←→  XINT(1+ni+MAXP1*nj, 1+nk+MAXP*nl)
!>   - loop over (i,j)           ←→  DO 7700 I / DO 7600 J
!>   - loop over (k,l)           ←→  DO 7500 K / DO 7400 L
!>   - lx formula                ←→  SOL(1) = (YINTI*ZINTJ - YINTJ*ZINTI)*XINT
!>   - dij_factor                ←→  FACI*CONJ(J)*CONK(K)*PNRM(K)*CONL(L)*PNRM(L)
!>                                   (EXPE is already inside gfull via F00)
!>
!> @param[in]    cpij    shell pair for electron 1 (bra=mu, ket=nu)
!> @param[in]    idij    primitive index within cpij
!> @param[in]    cpkl    shell pair for electron 2 (bra=lambda, ket=sigma)
!> @param[in]    idkl    primitive index within cpkl
!> @param[inout] socblk  accumulated block, shape (inao*jnao, 3)

 subroutine comp_soc_int2_prim(cpij, idij, cpkl, idkl, socblk)
!dir$ attributes inline :: comp_soc_int2_prim
    use ISO_FORTRAN_ENV, only: real64
    use mod_shell_tools, only: shpair_t
    use rys,             only: rys_root_t
    use constants,       only: CART_X, CART_Y, CART_Z, MAX_ANG => BAS_MXANG
    implicit none

    type(shpair_t),   intent(in)              :: cpij, cpkl
    integer,          intent(in)              :: idij, idkl
    real(real64), contiguous, intent(inout)   :: socblk(:, :)

    type(rys_root_t) :: ryscomp
    integer :: nroots_2e

    ! gfull(nj, ni, nl, nk, xyz, t) after VRR+HRR+unpack
    real(real64) :: gfull(0:MAX_ANG+1, 0:MAX_ANG+1, 0:MAX_ANG, 0:MAX_ANG, 3, (2*MAX_ANG+1)/2+2)

    ! derivatives for electron 1 (same arrays as comp_soc_int1_prim)
    integer,      parameter :: MAX_ANG_PAD = MAX_ANG + 1
    integer,      parameter :: MAX_NROOTS  = (2*MAX_ANG+1)/2 + 2
    real(real64) :: xyzin(0:MAX_ANG+1, 0:MAX_ANG+1, 3, MAX_NROOTS)
    real(real64) :: di(0:MAX_ANG_PAD, 0:MAX_ANG, 3, MAX_NROOTS)
    real(real64) :: dj(0:MAX_ANG_PAD, 0:MAX_ANG, 3, MAX_NROOTS)
!dir$ assume_aligned gfull  : 64
!dir$ assume_aligned socblk : 64

    ! loop indices
    integer  :: i, j, ij, jmax
    integer  :: k, l, lmax
    integer  :: nxi, nyi, nzi   ! Cartesian powers of function i (bra electron 1)
    integer  :: nxj, nyj, nzj   ! Cartesian powers of function j (ket electron 1)
    integer  :: nxk, nyk, nzk   ! Cartesian powers of function k (bra electron 2)
    integer  :: nxl, nyl, nzl   ! Cartesian powers of function l (ket electron 2)

    ! prefactor and SOC components
    real(real64) :: dij_fac    ! EXPE is inside gfull; this carries contraction coeff only
    real(real64) :: lx, ly, lz

    associate ( ppij => cpij%p(idij), ppkl => cpkl%p(idkl), &
                iang => cpij%iang,    jang => cpij%jang,     &
                inao => cpij%inao,    jnao => cpij%jnao,     &
                kang => cpkl%iang,    lang => cpkl%jang,     &
                knao => cpkl%inao,    lnao => cpkl%jnao )

    ! Number of Rys roots: total angular momentum of all four shells + 1
    ! GAMESS: NROOTS = MAXNM/2 + 1 where MAXNM = ILAM+JLAM+KLAM+LLAM+1
    nroots_2e      = (iang + jang + kang + lang + 1)/2 + 1
    ryscomp%nroots = nroots_2e

    ! Build the 4-index 1D integral table via Rys quadrature + VRR + HRR
    ! GAMESS equivalent: XYZ2E → XINT(N,M), YINT(N,M), ZINT(N,M)
    call QGaussRys2e(ryscomp, cpij, idij, cpkl, idkl, gfull)


!    if (iang==1 .and. jang==0 .and. kang==0 .and. idij==1 .and. idkl==1) then
!      write(6,'(a,3e20.12)') 'OQP gfull 00/10/01:', &
!      gfull(0,0,0,0,3,1), &   ! ZINT(1,1) = F00
!      gfull(0,1,0,0,1,1), &   ! XINT(2,1) = XC00
!      gfull(0,0,0,1,1,1)      ! XINT(1,2) = XCP00
!    endif
    ! Contraction prefactor for electron 1 primitive pair
    ! GAMESS: FACI*CONJ(J) are absorbed here; CONK*PNRM(K)*CONL*PNRM(L) go in (k,l) loop
    ! In OQP: pp%expfac already contains exp(-ai*aj/aa * |A-B|^2) * (pi/aa)^1.5
    !         EXPE is inside gfull (built into F00 = EXPE * w_t)
    !         So dij_fac here = 1.0 — all factors are already in gfull
    !         This mirrors how comp_soc_int1_prim uses dij = pp%expfac * TWOPI * pp%aa1
    ! EXPE is already inside gfull via F00 = EXPE*w(t) in QGaussRys2e.
    ! Contraction coefficients are applied in the outer loop (compute_som2e_ao).
    dij_fac = cpij%p(idij)%expfac * cpkl%p(idkl)%expfac! 1.0_real64

    ! --- Loop over subshell indices of electron 1 (GAMESS: DO 7700 I / DO 7600 J) ---
    ij   = 0
    jmax = jnao
    do i = 1, inao
        nxi = CART_X(i, iang);  nyi = CART_Y(i, iang);  nzi = CART_Z(i, iang)
        ! GAMESS: if IIEQJJ then JJMAX = I-1 (handled by cp%iandj in OQP)
        if (cpij%iandj) jmax = i - 1

        do j = 1, jmax
            nxj = CART_X(j, jang);  nyj = CART_Y(j, jang);  nzj = CART_Z(j, jang)
            ij = ij + 1

!            write(*,'(a,4i4)') 'DBG i j ij jmax=', i, j, ij, jmax
            ! Extract xyzin slice for this (i,j) pair from gfull:
            ! xyzin(nxl, nxk, 1, t) = gfull(nxj, nxi, nxl, nxk, 1, t)
            ! This replaces the XINT(NXX, MX) lookup in GAMESS

            ! --- Loop over subshell indices of electron 2 (GAMESS: DO 7500 K / DO 7400 L) ---
            lx = 0.0_real64;  ly = 0.0_real64;  lz = 0.0_real64
!       write(*,'(a,3i4,e14.6)') 'DBG i j ij lx=', i, j, ij, lx
            do k = 1, knao
                nxk = CART_X(k, kang);  nyk = CART_Y(k, kang);  nzk = CART_Z(k, kang)
                ! GAMESS: if KKEQLL then LLMAX = K
                lmax = lnao
                if (cpkl%iandj) lmax = k

                do l = 1, lmax
                    nxl = CART_X(l, lang);  nyl = CART_Y(l, lang);  nzl = CART_Z(l, lang)

                    ! Build xyzin for this (k,l) pair by taking the appropriate slice of gfull
                    ! GAMESS: MX=1+NX(K)+MAXP*NX(L), then XINT(NXX,MX) = gfull(nyj,nyi,nyl,nyk,2,t)
                    !
                    ! xyzin(nxl, nxk, xyz, t) = gfull(nxj, nxi, nxl, nxk, xyz, t)
                    ! but we need to build derivatives di, dj from xyzin first.
                    ! Here we directly use gfull elements in the Lx formula,
                    ! following GAMESS: SOL(1) = (YINTI(NYY,MY)*ZINTJ(NZZ,MZ)
                    !                           - YINTJ(NYY,MY)*ZINTI(NZZ,MZ)) * XINT(NXX,MX)
                    !
                    ! In OQP notation (after soc_xyz_ij on the (nxl,nxk) slice):
                    !   XINT(NXX,MX) = gfull(nxj, nxi, nxl, nxk, 1, t)   (no derivative)
                    !   YINTI(NYY,MY)= derivative of gfull w.r.t. bra-y on electron 1
                    !   ZINTJ(NZZ,MZ)= derivative w.r.t. ket-z on electron 1

                    ! Build xyzin for this (k,l) pair: three separate slices of gfull,
                    ! one per Cartesian component. Each component uses its own (nl,nk) index.
                    ! GAMESS: MX=1+NX(K)+MAXP*NX(L) selects col in XINT/YINT/ZINT.
                    !   xyzin(:,:,1,:) <- gfull(:,:, nxl, nxk, 1, :)
                    !   xyzin(:,:,2,:) <- gfull(:,:, nyl, nyk, 2, :)
                    !   xyzin(:,:,3,:) <- gfull(:,:, nzl, nzk, 3, :)
                    xyzin(0:jang+1, 0:iang+1, 1, 1:nroots_2e) = &
                        gfull(0:jang+1, 0:iang+1, nxl, nxk, 1, 1:nroots_2e)
                    xyzin(0:jang+1, 0:iang+1, 2, 1:nroots_2e) = &
                        gfull(0:jang+1, 0:iang+1, nyl, nyk, 2, 1:nroots_2e)
                    xyzin(0:jang+1, 0:iang+1, 3, 1:nroots_2e) = &
                        gfull(0:jang+1, 0:iang+1, nzl, nzk, 3, 1:nroots_2e)

                    ! Build derivatives di(m,n) = n*xyzin(m,n-1) - 2*ai*xyzin(m,n+1)  [bra]
                    !              and dj(m,n) = m*xyzin(m-1,n) - 2*aj*xyzin(m+1,n)  [ket]
                    ! GAMESS: YINTI(NYY,MY) = NI*YINT(NI-1,MY) - 2*AI*YINT(NI+1,MY)
                    call soc_xyz_ij(xyzin, iang, jang, ppij%ai, ppij%aj, nroots_2e, di, dj)

!if (iang==1 .and. jang==0 .and. kang==0 .and. idij==1 .and. idkl==1) then
!  write(6,'(a,4i3,6e16.8)') 'OQP di/dj i/j/k/l=', i,j,k,l, &
!    di(nyj,nyi,2,1), dj(nyj,nyi,2,1), &
!    di(nzj,nzi,3,1), dj(nzj,nzi,3,1), &
!    di(nxj,nxi,1,1), dj(nxj,nxi,1,1)
!endif

                    ! --- Lx, Ly, Lz (GAMESS: SOL(1), SOL(2), SOL(3)) ---
                    ! SOL(1) = (YINTI(NYY,MY)*ZINTJ(NZZ,MZ) - YINTJ(NYY,MY)*ZINTI(NZZ,MZ))
                    !        * XINT(NXX,MX)
                    ! In OQP:
                    !   XINT(NXX,MX)   = sum_t xyzin(nxj,nxi,1,t)   [no derivative]
                    !   YINTI(NYY,MY)  = sum_t di(nyj,nyi,2,t)
                    !   ZINTJ(NZZ,MZ)  = sum_t dj(nzj,nzi,3,t)
                    !   YINTJ(NYY,MY)  = sum_t dj(nyj,nyi,2,t)
                    !   ZINTI(NZZ,MZ)  = sum_t di(nzj,nzi,3,t)

                    lx = lx + sum( xyzin(nxj, nxi, 1, 1:nroots_2e) &
                                 *    di(nyj, nyi, 2, 1:nroots_2e)  &
                                 *    dj(nzj, nzi, 3, 1:nroots_2e) )&
                             - sum( xyzin(nxj, nxi, 1, 1:nroots_2e) &
                                 *    di(nzj, nzi, 3, 1:nroots_2e)  &
                                 *    dj(nyj, nyi, 2, 1:nroots_2e) )

                    ly = ly + sum( xyzin(nyj, nyi, 2, 1:nroots_2e) &
                                 *    di(nzj, nzi, 3, 1:nroots_2e)  &
                                 *    dj(nxj, nxi, 1, 1:nroots_2e) )&
                             - sum( xyzin(nyj, nyi, 2, 1:nroots_2e) &
                                 *    di(nxj, nxi, 1, 1:nroots_2e)  &
                                 *    dj(nzj, nzi, 3, 1:nroots_2e) )

                    lz = lz + sum( xyzin(nzj, nzi, 3, 1:nroots_2e) &
                                 *    di(nxj, nxi, 1, 1:nroots_2e)  &
                                 *    dj(nyj, nyi, 2, 1:nroots_2e) )&
                             - sum( xyzin(nzj, nzi, 3, 1:nroots_2e) &
                                 *    di(nyj, nyi, 2, 1:nroots_2e)  &
                                 *    dj(nxj, nxi, 1, 1:nroots_2e) )


!                    if (iang==1 .and. jang==0 .and. kang==0 .and. idij==1 .and. idkl==1) then
!                      write(6,'(a,4i3,3e20.12)') 'OQP lx/ly/lz i/j/k/l=', i,j,k,l, lx,ly,lz
!                    endif
                end do  ! l
            end do  ! k
!if (idij==1 .and. idkl==1 .and. i==2 .and. j==1) then
!    write(*,'(a,3e14.6)') 'OQP lx ly lz:', lx, ly, lz
!end if
            ! Accumulate into socblk
            ! GAMESS: SO2AO -= TDENFC * SOL  where TDENFC = FACK*CONL*PNRM(L)
            ! In OQP: dij_fac carries the primitive prefactor; contraction
            !         weights from p_kl are applied in the outer loop (compute_som2e_ao)
            socblk(ij, 1) = socblk(ij, 1) - dij_fac * lx
            socblk(ij, 2) = socblk(ij, 2) - dij_fac * ly
            socblk(ij, 3) = socblk(ij, 3) - dij_fac * lz

        end do  ! j
    end do  ! i

    end associate

 end subroutine comp_soc_int2_prim










! 2e part ends here. 




END MODULE
