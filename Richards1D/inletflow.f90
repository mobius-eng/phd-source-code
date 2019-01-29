MODULE INLETFLOW_MODULE
    
    USE GENERAL
    USE GROWING_ARRAY_MOD
    IMPLICIT NONE
    
    ! INLETFLOW
    !   Base type for all inlet flow types. Defaults to no flow.
    !
    TYPE INLETFLOW
    CONTAINS
        PROCEDURE :: Q => DEFAULT_Q    
    END TYPE
    
    
    ! CONST_INLETFLOW
    !   Constant inlet flow
    !
    TYPE, EXTENDS(INLETFLOW) :: CONST_INLETFLOW
        REAL(DP) :: CONSTQ
    CONTAINS
        PROCEDURE :: Q => CONST_INLETFLOW_Q
    END TYPE
    
    
    ! INTERMITTENT_INLETFLOW
    !   Intermittent flow with flow on for TON period at QON flow rate
    !   and flow off for TOFF time period.
    !
    TYPE, EXTENDS(INLETFLOW) :: INTERMITTENT_INLETFLOW
        REAL(DP) :: QON, TON, TOFF
    CONTAINS
        PROCEDURE :: Q => INTERMITTENT_INLETFLOW_Q
    END TYPE
    
    TYPE, EXTENDS(INLETFLOW) :: NOISY_INLETFLOW
        REAL(DP) :: QBASE, NOISE, PERIOD
        TYPE(GROWING_ARRAY) :: QHIST
    CONTAINS
        PROCEDURE :: Q => NOISY_INLETFLOW_Q
    END TYPE
    
    

    CONTAINS
    
    REAL(DP) FUNCTION DEFAULT_Q(QIN, T)
    ! DEFAULT_Q
    !   Default flow function Q for INLETFLOW: no flow
    !   
    CLASS(INLETFLOW), INTENT(INOUT) :: QIN
    REAL(DP), INTENT(IN) :: T
    DEFAULT_Q = 0D0
    END FUNCTION DEFAULT_Q
    
    
    REAL(DP) FUNCTION CONST_INLETFLOW_Q(QIN, T)
    ! CONST_INLETFLOW
    !   Inlet flow function for the flow at constant rate
    !
    CLASS(CONST_INLETFLOW), INTENT(INOUT) :: QIN
    REAL(DP), INTENT(IN) :: T
    CONST_INLETFLOW_Q = QIN % CONSTQ
    END FUNCTION CONST_INLETFLOW_Q
    
    
    REAL(DP) FUNCTION INTERMITTENT_INLETFLOW_Q(QIN, T) RESULT(Q)
    ! INTERMITTENT_INLETFLOW
    !   "Square wave" fow function
    CLASS(INTERMITTENT_INLETFLOW), INTENT(INOUT) :: QIN
    REAL(DP), INTENT(IN) :: T
    
    REAL(DP) :: TPERIOD, TRES
    ! Reduce T to the single period: 0 <= TRES < TON+TOFF
    TPERIOD = QIN % TON + QIN % TOFF
    TRES = MOD(T, TPERIOD)
    ! Find where in the period
    IF (TRES > QIN % TON) THEN
        Q = 0D0
    ELSE
        Q = QIN % QON
    END IF
    END FUNCTION INTERMITTENT_INLETFLOW_Q
    
    
    FUNCTION MK_NOISY_INLETFLOW(QBASE, NOISE, PERIOD, SZ) RESULT(QN)
        REAL(DP), INTENT(IN) :: QBASE, NOISE, PERIOD
        INTEGER, INTENT(IN), OPTIONAL :: SZ
        TYPE(NOISY_INLETFLOW) :: QN
        
        INTEGER :: ISZ
        
        IF (.NOT. PRESENT(SZ)) THEN
            ISZ = 1000
        ELSE
            ISZ = SZ
        END IF
        
        QN%QBASE = QBASE
        QN%NOISE = NOISE
        QN%PERIOD = PERIOD
        QN%QHIST = MKGRARRAY(ISZ, I0 = 0)
        CALL RANDOM_SEED
    END FUNCTION MK_NOISY_INLETFLOW
    
    REAL(DP) FUNCTION NOISY_INLETFLOW_Q(QIN, T) RESULT(Q)
        CLASS(NOISY_INLETFLOW), INTENT(INOUT) :: QIN
        REAL(DP), INTENT(IN) :: T
        
        INTEGER :: PERIODNUMBER
        REAL(DP) :: QHIST, RND
        LOGICAL :: FOUND
        
        PERIODNUMBER = INT(T / QIN%PERIOD)
        QHIST = QIN%QHIST%GET(PERIODNUMBER, FOUND)
        IF (FOUND) THEN
            Q = QHIST
            RETURN
        ELSE
            CALL RANDOM_NUMBER(RND)
            RND = (RND - 0.5D0) * 2D0
            Q = QIN%QBASE * (1D0 + QIN%NOISE * RND)
            CALL QIN%QHIST%SET(PERIODNUMBER, Q)
            RETURN
        END IF
    END FUNCTION NOISY_INLETFLOW_Q
    
END MODULE INLETFLOW_MODULE