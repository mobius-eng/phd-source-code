MODULE CSTR

USE ENV

IMPLICIT NONE

TYPE CSTR_STATE
    REAL(WREAL) :: C, X
END TYPE CSTR_STATE


TYPE CSTR_CONFIG
    REAL(WREAL) :: CFLOW, CREACT
    PROCEDURE(INLET_FUNCTION), NOPASS, POINTER :: C_INLET
    PROCEDURE(SOURCE_TERM_PROCEDURE), NOPASS, POINTER :: SOURCE_TERM
END TYPE CSTR_CONFIG


ABSTRACT INTERFACE

    FUNCTION INLET_FUNCTION(TIME)
    
    IMPORT :: WREAL
    REAL(WREAL), INTENT(IN) :: TIME
    REAL(WREAL) :: INLET_FUNCTION
    
    END FUNCTION INLET_FUNCTION
    
    
    SUBROUTINE SOURCE_TERM_PROCEDURE(STATE, F, FC, FX)
    
    IMPORT :: WREAL, CSTR_STATE
    
    TYPE(CSTR_STATE), INTENT(IN) :: STATE
    REAL(WREAL), INTENT(OUT) :: F, FC, FX
    
    END SUBROUTINE SOURCE_TERM_PROCEDURE
    
END INTERFACE



CONTAINS


PURE FUNCTION CSTR_NEWTON_STEP(A, B) RESULT (Y)

REAL(WREAL), DIMENSION(1:2, 1:2), INTENT(IN) :: A
REAL(WREAL), DIMENSION(1:2), INTENT(IN) :: B
TYPE(CSTR_STATE) :: Y

REAL(WREAL) :: DETA, DETA1, DETA2

DETA  = A(1,1) * A(2,2) - A(1,2) * A(2,1)
DETA1 = B(1) * A(2,2) - B(2) * A(1,2);
DETA2 = A(1,1) * B(2) - B(1) * A(2,1);

Y%C = DETA1 / DETA
Y%X = DETA2 / DETA

END FUNCTION CSTR_NEWTON_STEP


SUBROUTINE FORM_COEFF_FOR_NEWTON_STEP(CONF, DT, B0, STATE, A, B)

TYPE(CSTR_CONFIG), INTENT(IN) :: CONF
REAL(WREAL), INTENT(IN) :: DT
TYPE(CSTR_STATE), INTENT(IN) :: STATE
REAL(WREAL), DIMENSION(1:2,1:2), INTENT(OUT) :: A
REAL(WREAL), DIMENSION(1:2), INTENT(OUT) :: B
REAL(WREAL), DIMENSION(1:2), INTENT(IN) :: B0

REAL(WREAL) :: F, FC, FX, PHI

CALL CONF%SOURCE_TERM(STATE, F, FC, FX)
PHI = F - FC * STATE%C - FX * STATE%X
B(1) = B0(1) - DT * PHI / 2
B(2) = B0(2) + CONF%CREACT * DT * PHI / 2
A(1,1) = 1.0_WREAL + CONF%CFLOW * DT / 2 + DT * FC / 2
A(1,2) = DT * FX / 2;
A(2,1) = - CONF%CREACT * DT * FC / 2
A(2,2) = 1.0_WREAL - CONF%CREACT * DT * FX / 2

END SUBROUTINE FORM_COEFF_FOR_NEWTON_STEP


FUNCTION CSTR_TIME_STEP(CONF, STATE, TIME, DT, ATOL, RTOL, MAXITER) RESULT (NEWSTATE)

TYPE(CSTR_CONFIG), INTENT(IN) :: CONF
TYPE(CSTR_STATE), INTENT(IN) :: STATE
REAL(WREAL), INTENT(IN) :: TIME, DT, ATOL, RTOL
INTEGER(WINT), INTENT(IN) :: MAXITER

REAL(WREAL) :: F, FC, FX, CIN0, CIN1
REAL(WREAL), DIMENSION(1:2,1:2) :: A
REAL(WREAL), DIMENSION(1:2) :: B, B0
INTEGER(WINT) :: ITER
TYPE(CSTR_STATE) :: OLDSTATE, NEWSTATE

CALL CONF%SOURCE_TERM(STATE, F, FC, FX)

CIN0 = CONF%C_INLET(TIME)
CIN1 = CONF%C_INLET(TIME + DT)

B0(1) = CONF%CFLOW * DT * (CIN0 + CIN1) / 2 + (1 - CONF%CFLOW * DT / 2) * STATE%C - DT * F / 2
B0(2) = STATE%X + CONF%CREACT * DT * F
NEWSTATE = STATE

DO ITER = 1, MAXITER
    OLDSTATE = NEWSTATE
    CALL FORM_COEFF_FOR_NEWTON_STEP(CONF, DT, B0, NEWSTATE, A, B)
    NEWSTATE = CSTR_NEWTON_STEP(A, B)
    IF (ABS(NEWSTATE%C - OLDSTATE%C) < MAX(ATOL, RTOL * ABS(NEWSTATE%C))) RETURN
END DO

NEWSTATE%C = -1.0_WREAL
NEWSTATE%X = -1.0_WREAL

END FUNCTION CSTR_TIME_STEP

SUBROUTINE RUN_CSTR_SIMULATION(CONF, TIME0, DT, STATES, ATOL, RTOL, MAXITER)

TYPE(CSTR_CONFIG), INTENT(IN) :: CONF
REAL(WREAL), INTENT(IN) :: TIME0, DT, ATOL, RTOL
TYPE(CSTR_STATE), DIMENSION(0:), INTENT(INOUT) :: STATES
INTEGER(WINT), INTENT(IN) :: MAXITER

INTEGER(WINT) :: ISTEP
REAL(WREAL) :: TIME

TIME = TIME0

DO ISTEP = 1, SIZE(STATES, 1)
    STATES(ISTEP) = CSTR_TIME_STEP(CONF, STATES(ISTEP-1), TIME, DT, ATOL, RTOL, MAXITER)
    TIME = TIME + DT
END DO

END SUBROUTINE RUN_CSTR_SIMULATION

END MODULE CSTR