MODULE GROWING_ARRAY_MOD
    ! USE GENERAL
    
    IMPLICIT NONE
    INTEGER, PARAMETER, PRIVATE :: DP = KIND(1D0)
    
    TYPE GROWING_ARRAY
        INTEGER :: SIZE, UBOUND, LBOUND
        REAL(DP), DIMENSION(:,:), ALLOCATABLE :: BUFFER
        TYPE(GROWING_ARRAY), POINTER :: NEXT => NULL()
    CONTAINS
        PROCEDURE, PASS(THIS) :: GET => GRARR_GET
        PROCEDURE, PASS(THIS) :: SET => GRARR_SET
    END TYPE
    
    CONTAINS
    
    FUNCTION MKGRARRAY(CHSIZE, ISIZE, I0) RESULT (ARR)
        INTEGER, INTENT(IN) :: CHSIZE
        INTEGER, INTENT(IN), OPTIONAL :: ISIZE, I0
        TYPE(GROWING_ARRAY),  TARGET :: ARR
        
        INTEGER :: ICUR, LBOUNDMIN, UBOUNDMAX
        TYPE(GROWING_ARRAY), POINTER :: CURARR
        
        IF (.NOT. PRESENT(I0)) THEN
            LBOUNDMIN = 1
        ELSE
            LBOUNDMIN = I0
        END IF
        
        IF (.NOT. PRESENT(ISIZE)) THEN
            UBOUNDMAX = CHSIZE + LBOUNDMIN - 1
        ELSE
            UBOUNDMAX = ISIZE + LBOUNDMIN - 1
        END IF
        
        CURARR => ARR
        ! Place it before the first element
        ICUR = LBOUNDMIN - 1
        
        DO WHILE (ICUR + 1 < UBOUNDMAX)
            CURARR%SIZE = CHSIZE
            CURARR%LBOUND = ICUR + 1
            CURARR%UBOUND = CHSIZE + ICUR
            ALLOCATE (CURARR%BUFFER(2, CURARR%LBOUND : CURARR%UBOUND), SOURCE=0D0)
            ! PLACE IT BEFORE THE NEXT ELEMENT
            ICUR = ICUR + CHSIZE
            ! Check if need another chunck
            IF (ICUR + 1 < UBOUNDMAX) THEN
                ALLOCATE (CURARR%NEXT)
                CURARR => CURARR%NEXT
            END IF
        END DO
        
    END FUNCTION MKGRARRAY
    
    FUNCTION GRARR_GET(THIS, I, PRESENT) RESULT(EL)
        CLASS(GROWING_ARRAY), INTENT(IN), TARGET :: THIS
        INTEGER, INTENT(IN) :: I
        LOGICAL, INTENT(OUT) :: PRESENT
        REAL(DP) :: EL, ISELSET
        
        TYPE(GROWING_ARRAY), POINTER :: PARR
        
        PARR => THIS
        
        EL = 0D0
        PRESENT = .FALSE.
        DO
            IF(I <= PARR%UBOUND) THEN
                EL = PARR%BUFFER(1, I)
                ISELSET = PARR%BUFFER(2,I)
                IF (ISELSET > 0D0) PRESENT = .TRUE.
                RETURN
            END IF
            PARR => PARR%NEXT
            IF (.NOT. ASSOCIATED(PARR)) RETURN
        END DO
    END FUNCTION
    
    SUBROUTINE GRARR_SET(THIS, I, V)
        CLASS(GROWING_ARRAY), INTENT(INOUT), TARGET :: THIS
        INTEGER, INTENT(IN) :: I
        REAL(DP), INTENT(IN) :: V
        
        TYPE(GROWING_ARRAY), POINTER :: PARR
        
        PARR => THIS
        
        DO
            IF (I <= PARR%UBOUND) THEN
                PARR%BUFFER(1,I) = V
                PARR%BUFFER(2,I) = 1D0
                RETURN
            END IF
            IF (.NOT. ASSOCIATED(PARR%NEXT)) THEN
                ALLOCATE (PARR%NEXT)
                PARR%NEXT%SIZE = PARR%SIZE
                PARR%NEXT%LBOUND = PARR%UBOUND+1
                PARR%NEXT%UBOUND = PARR%UBOUND + PARR%NEXT%SIZE
                ALLOCATE (PARR%NEXT%BUFFER(2, PARR%NEXT%LBOUND : PARR%NEXT%UBOUND), SOURCE=0D0)
            END IF
            PARR => PARR%NEXT
        END DO
    END SUBROUTINE GRARR_SET
    
    RECURSIVE SUBROUTINE DEL_GROWING_ARRAY(GRARR)
        TYPE(GROWING_ARRAY), INTENT(INOUT) :: GRARR
        
        IF (ASSOCIATED(GRARR%NEXT)) THEN
            CALL DEL_GROWING_ARRAY(GRARR%NEXT)
        END IF
        
        DEALLOCATE (GRARR%BUFFER)
        GRARR%NEXT => NULL()
        GRARR%SIZE = 0
        GRARR%UBOUND = 0
        GRARR%LBOUND = 0
        
    END SUBROUTINE
    
    
END MODULE GROWING_ARRAY_MOD