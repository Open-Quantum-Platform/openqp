MODULE mod_shell_tools

 USE precision, only: dp
 USE basis_tools, ONLY: basis_set

 TYPE shell_t
     !INTEGER :: shid, atid, ig1, ig2, ang, minbf, maxbf, locao, nao
     INTEGER :: shid       !< shell ID in a basis set
     INTEGER :: atid       !< ID of the corresponding atom
     INTEGER :: ig1        !< index of the first primitive of the shell in total basis set
     INTEGER :: ig2        !< index of the last primitive of the shell in total basis set
     INTEGER :: ang        !< angular momentum + 1
     INTEGER :: locao      !< position of shell's first basis function in the basis set
     INTEGER :: nao        !< number of basis functions (atomic orbitals) in shell
     REAL(kind=dp) :: r(3) !< shell origin
     CONTAINS
     PROCEDURE :: fetch_by_id => bas_set_indices
 END TYPE

!< @brief Data structure to hold data of primitive shell pair
 TYPE primpair_t
    REAL(kind=dp) :: r(3)      !< center of charge coordinates
    REAL(kind=dp) :: aa        !< sum of exponents
    REAL(kind=dp) :: aa1       !< reverse sum of exponents
    REAL(kind=dp) :: ai        !< exponent of the first shell
    REAL(kind=dp) :: aj        !< exponent of the second shell
    REAL(kind=dp) :: expfac    !< exponential prefactor
 END TYPE

!< @brief Data structure to hold data of contracted shell pair
!   I think AoS for primitive pairs would be OK as the only computationally
!   intensive part involving primitive pairs formation is not the limiting step.
!   There is small benefit in vectorizing 1e integrals: most of them are in fact
!   memory bandwidth limited, so it is better to focus on smaller CPU cache footprint.
 TYPE shpair_t
    REAL(kind=dp) :: ri(3)              !< coordinates of the first shell
    REAL(kind=dp) :: pad1
    REAL(kind=dp) :: rj(3)              !< coordinates of the second shell
    REAL(kind=dp) :: pad2
    INTEGER :: iang                     !< angular momentum of the first shell
    INTEGER :: jang                     !< angular momentum of the second shell
    INTEGER :: inao                     !< number of B.F. in the first shell
    INTEGER :: jnao                     !< number of B.F. in the second shell
    INTEGER :: isder                    !< `isder/=0` if derivatives are computed
    INTEGER :: numpairs                 !< number of primitives in shell pair block
    INTEGER :: nroots
    LOGICAL :: iandj                    !< true if first shell is same to second one
    TYPE(primpair_t), ALLOCATABLE :: p(:) !< array of primitive pairs data
!dir$ attributes align : 64 :: p
    CONTAINS
    PROCEDURE :: alloc => shell_pair_alloc
    PROCEDURE :: alloc2 => shell_pair_alloc2
    PROCEDURE :: shell_pair
    PROCEDURE :: shell_pair2
 END TYPE

 PRIVATE

 PUBLIC shell_t
 PUBLIC shpair_t

! PUBLIC shell_pair
! PUBLIC shell_pair_alloc
! PUBLIC bas_set_indices

CONTAINS

!> @brief Copy shell info from basis set to shell_t variable
!> @param[in]  shid     index of shell in a basis set
!> @param[out] shinfo   shell data
!
!> @author   Vladimir Mironov
!
!     REVISION HISTORY:
!> @date _Sep, 2018_ Initial release
!
 SUBROUTINE bas_set_indices(shinfo, basis, shid)
!dir$ attributes inline :: bas_set_indices
    INTEGER, INTENT(IN) :: shid
    type(basis_set), INTENT(IN) :: basis
    CLASS(shell_t), INTENT(INOUT) :: shinfo

    shinfo%shid  = shid
    shinfo%atid  = basis%origin(shid)
    shinfo%ig1   = basis%g_offset(shid)
    shinfo%ig2   = basis%g_offset(shid)+basis%ncontr(shid)-1
    shinfo%ang   = basis%am(shid)
    shinfo%locao = basis%ao_offset(shid)
    shinfo%nao   = basis%naos(shid)
    shinfo%r     = basis%atoms%xyz(:,shinfo%atid)

 END SUBROUTINE bas_set_indices

!--------------------------------------------------------------------------------

!> @brief Allocate shell pair array of primitives
!> @param[out] sp     shell pair data
!> @param[in]  basis  basis set
!
!> @author   Vladimir Mironov
!
!     REVISION HISTORY:
!> @date _Oct, 2018_ Initial release
!
 SUBROUTINE shell_pair_alloc(sp, basis)
    CLASS(shpair_t), INTENT(INOUT) :: sp
    type(basis_set), INTENT(IN) :: basis

    IF(allocated(sp%p).AND.ubound(sp%p,1)<basis%mxcontr**2) DEALLOCATE(sp%p)

    IF (.NOT.allocated(sp%p)) ALLOCATE(sp%p(basis%mxcontr**2))

 END SUBROUTINE

!--------------------------------------------------------------------------------

!> @brief Allocate shell pair array of primitives
!> @param[out] sp     shell pair data
!> @param[in]  basis1 first basis set
!> @param[in]  basis2 second basis set
!
!> @author   Igor S. Gerasimov
!
!     REVISION HISTORY:
!> @date _Oct, 2022_ Initial release
!
 SUBROUTINE shell_pair_alloc2(sp, basis1, basis2)
    CLASS(shpair_t), INTENT(INOUT) :: sp
    type(basis_set), INTENT(IN) :: basis1, basis2

    IF(allocated(sp%p).AND.ubound(sp%p,1)<basis1%mxcontr*basis2%mxcontr) DEALLOCATE(sp%p)

    IF (.NOT.allocated(sp%p)) ALLOCATE(sp%p(basis1%mxcontr*basis2%mxcontr))

 END SUBROUTINE

!--------------------------------------------------------------------------------

!> @brief Generate data for pairwise electronic distributions
!> @param[in]  shi      index of the first shell in a basis set
!> @param[in]  shj      index of the second shell in a basis set
!> @param[in]  tol      cut-off for integrals
!> @param[out] pair     shell pair data
!> @param[in]  dup      [optional] if `.true.` - multiply prefactors by 2 when `shi==shj`
!
!> @author   Vladimir Mironov
!
!     REVISION HISTORY:
!> @date _Oct, 2018_ Initial release
!
 SUBROUTINE shell_pair(pair, basis, shi, shj, tol, dup)
!dir$ attributes inline :: shell_pair
    CLASS(shpair_t), INTENT(INOUT) :: pair
    type(basis_set), INTENT(IN) :: basis
    TYPE(shell_t), INTENT(IN) :: shi, shj
    REAL(kind=dp), INTENT(IN) :: tol
    LOGICAL, INTENT(IN), OPTIONAL :: dup

    REAL(kind=dp) :: rr, cci, ccj, aa, fac, ai, aj
    INTEGER :: ig, jg, jgmax, ij
    LOGICAL :: iandj, duplicate
    duplicate = .TRUE.
    if (present(dup)) duplicate = dup

    rr = sum( (shi%r-shj%r)**2 )

    pair%ri = shi%r
    pair%rj = shj%r

    pair%nroots = (shi%ang+shj%ang+1)/2 + 1

    iandj = shi%shid==shj%shid
    pair%iandj = iandj
    iandj = iandj.AND.duplicate

    pair%iang = shi%ang
    pair%jang = shj%ang

    pair%inao = shi%nao
    pair%jnao = shj%nao

!   I primitive
    ij = 0
    jgmax = shj%ig2
    DO ig = shi%ig1, shi%ig2

        ai = basis%ex(ig)
        cci = basis%cc(ig)

!       J primitive
        IF (iandj) jgmax = ig
        DO jg = shj%ig1, jgmax

            aj = basis%ex(jg)
            ccj = basis%cc(jg)

            aa = ai+aj

            IF (ai*aj*rr > tol*aa) CYCLE

            ij = ij + 1

            ASSOCIATE (pp => pair%p(ij))
                pp%ai = ai
                pp%aj = aj

                pp%aa = aa
                pp%aa1 = 1/aa

                pp%r = (ai*pair%ri + aj*pair%rj) * pp%aa1

                fac = cci*ccj*exp(-ai*aj*pp%aa1*rr)
                pp%expfac = fac
                IF (iandj) pp%expfac = 2*fac
            END ASSOCIATE

        END DO
        IF (iandj) pair%p(ij)%expfac = fac
    END DO

    pair%numpairs = ij

 END SUBROUTINE

!--------------------------------------------------------------------------------

!> @brief Generate data for pairwise electronic distributions for basis sets overlapping
!> @param[out] pair     shell pair data
!> @param[in]  shi      index of the first shell in a first basis set
!> @param[in]  shj      index of the second shell in a second basis set
!> @param[in]  tol      cut-off for integrals
!
!> @author   Igor S. Gerasimov
!
!     REVISION HISTORY:
!> @date _Oct, 2022_ Initial release
!
 SUBROUTINE shell_pair2(pair, basis1, basis2, shi, shj, tol)
!dir$ attributes inline :: shell_pair2
    CLASS(shpair_t), INTENT(INOUT) :: pair
    type(basis_set), INTENT(IN) :: basis1, basis2
    TYPE(shell_t), INTENT(IN) :: shi, shj
    REAL(kind=dp), INTENT(IN) :: tol

    REAL(kind=dp) :: rr, cci, ccj, aa, fac, ai, aj
    INTEGER :: ig, jg, ij

    rr = sum( (shi%r-shj%r)**2 )

    pair%ri = shi%r
    pair%rj = shj%r

    pair%nroots = (shi%ang+shj%ang+1)/2 + 1

    pair%iandj = .false.

    pair%iang = shi%ang
    pair%jang = shj%ang

    pair%inao = shi%nao
    pair%jnao = shj%nao

!   I primitive
    ij = 0
    DO ig = shi%ig1, shi%ig2

        ai = basis1%ex(ig)
        cci = basis1%cc(ig)

!       J primitive
        DO jg = shj%ig1, shj%ig2

            aj = basis2%ex(jg)
            ccj = basis2%cc(jg)

            aa = ai+aj

            IF (ai*aj*rr > tol*aa) CYCLE

            ij = ij + 1

            ASSOCIATE (pp => pair%p(ij))
                pp%ai = ai
                pp%aj = aj

                pp%aa = aa
                pp%aa1 = 1/aa

                pp%r = (ai*pair%ri + aj*pair%rj) * pp%aa1

                fac = cci*ccj*exp(-ai*aj*pp%aa1*rr)
                pp%expfac = fac
            END ASSOCIATE

        END DO
    END DO

    pair%numpairs = ij

 END SUBROUTINE

!--------------------------------------------------------------------------------
END MODULE
