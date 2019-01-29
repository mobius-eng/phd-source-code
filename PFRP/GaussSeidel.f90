MODULE GAUSS_SEIDEL

USE ENV

IMPLICIT NONE

PRIVATE

TYPE, PUBLIC :: GS_PROBLEM(N)
    INTEGER(WINT), LEN :: N
    REAL(WREAL), DIMENSION(0:N) :: AC2, AC1, ACX, ACM, AXC, AXM, BC, BX
    REAL(WREAL) :: TOL, ERR
    INTEGER(WINT) :: MAXITER, ITER
END TYPE GS_PROBLEM

INTERFACE SOLVE
    PROCEDURE :: SOLVE_1, SOLVE_2
END INTERFACE SOLVE

PUBLIC :: SOLVE, SOLVE_1, SOLVE_2

CONTAINS

SUBROUTINE STEP(C, X, ERR, AC2, AC1, ACX, ACM, AXC, AXM, BC, BX)

REAL(WREAL), DIMENSION(0:), INTENT(INOUT) :: C, X
REAL(WREAL), INTENT(OUT) :: ERR
REAL(WREAL), DIMENSION(0:), INTENT(IN) :: AC2, AC1, ACX, ACM, AXC, AXM, BC, BX

REAL(WREAL) :: TMP, TMP2
INTEGER(WINT) :: I

! Init error as zero. It will only increase later.
ERR = 0.0_WREAL
! Top boundary node: C = 1 is specified. Compute X
C(0) = 1.0_WREAL
X(0) = MIN( MAX( AXM(0) * (AXC(0) * C(0) + BX(0)), 0.0_WREAL ), 1.0_WREAL )
! First inner node: use linear derivative approximation in convective term
TMP = C(1)
C(1) = MAX( ACM(1) * (AC1(1)*C(0) + ACX(1) * X(1) + BC(1)), 0.0_WREAL )
X(1) = MIN( MAX( AXM(1) * (AXC(1) * C(1) + BX(1)), 0.0_WREAL ), 1.0_WREAL )
ERR = MAX( ERR, ABS(TMP - C(1)) )
! Other nodes: use quadratic derivative approximation
DO I = 2, SIZE(C)-1
    TMP = C(I)
    TMP2 = ACM(I) * (AC2(I) * C(I-2) + AC1(I) * C(I-1) + ACX(I) * X(I) + BC(I))
    C(I) = MAX( TMP2, 0.0_WREAL )
    X(I) = MIN( MAX( AXM(I) * (AXC(I) * C(I) + BX(I)), 0.0_WREAL ), 1.0_WREAL )
    ERR = MAX( ERR, ABS(TMP - C(I)) )
END DO
END SUBROUTINE STEP


SUBROUTINE SOLVE_1(C, X, AC2, AC1, ACX, ACM, AXC, AXM, BC, BX, TOL, MAXITER, ERR, ITER)

REAL(WREAL), DIMENSION(0:), INTENT(INOUT) :: C, X
REAL(WREAL), DIMENSION(0:), INTENT(IN) :: AC2, AC1, ACX, ACM, AXC, AXM, BC, BX
REAL(WREAL), INTENT(IN) :: TOL
INTEGER(WINT), INTENT(IN) :: MAXITER
REAL(WREAL), INTENT(OUT) :: ERR
INTEGER(WINT), INTENT(OUT) :: ITER

DO ITER = 1, MAXITER
    CALL STEP(C, X, ERR, AC2, AC1, ACX, ACM, AXC, AXM, BC, BX)
    IF (ERR < TOL) THEN
        RETURN
    END IF
END DO
ERR = -1.0
RETURN

END SUBROUTINE SOLVE_1


SUBROUTINE SOLVE_2(C, X, COEFFS)

TYPE(GS_PROBLEM(*)), INTENT(INOUT) :: COEFFS
REAL(WREAL), DIMENSION(0:COEFFS%N), INTENT(INOUT) :: C, X
INTEGER(WINT) :: ITER

CALL SOLVE_1(C, X, COEFFS%AC2, COEFFS%AC1, COEFFS%ACX, COEFFS%ACM, &
    COEFFS%AXC, COEFFS%AXM, COEFFS%BC, COEFFS%BX, COEFFS%TOL, COEFFS%MAXITER, COEFFS%ERR, ITER)

COEFFS%ITER = COEFFS%ITER + ITER

END SUBROUTINE


END MODULE GAUSS_SEIDEL