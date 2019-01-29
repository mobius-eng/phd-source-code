! The interface to the library

MODULE RUNRICHARDS

CONTAINS
    
    INTEGER FUNCTION RUN_RICHARDSUC(KS, L, A, N, M, DTHETA, QIN, ZMAX, NZ, TMAX, NT, S0, &
        TOUT, ZOUT, SATOUT, QINOUT, QOUTOUT) RESULT(ERR)
    !DEC$ ATTRIBUTES DLLEXPORT, ALIAS : 'RUN_RICHARDSUC' :: RUN_RICHARDSUC
    USE GENERAL
    USE INLETFLOW_MODULE
    USE UNSATURATED
    USE DARCY_MODULE
    IMPLICIT NONE
    
    REAL(DP), INTENT(IN) :: KS, L, A, N, M, DTHETA, QIN, ZMAX, TMAX, S0
    INTEGER, INTENT(IN) :: NZ, NT
    REAL(DP), INTENT(OUT) :: TOUT(NT), ZOUT(NZ), SATOUT(NZ, NT), QINOUT(NT), QOUTOUT(NT)
    
    TYPE(VANGENUCHTEN) :: UMDL
    TYPE(DARCY) :: DMDL
    TYPE(CONST_INLETFLOW) :: QFLOW
    REAL(DP), DIMENSION(:), ALLOCATABLE :: SAT0
    
    UMDL = MKVANGENUCHTEN(L, A, N, M)
    DMDL = MKDARCY(NZ, KS, DTHETA, UMDL)
    QFLOW = CONST_INLETFLOW(QIN)
    ERR = 0
    
    ALLOCATE (SAT0(NZ), STAT = ERR)
    IF (ERR > 0) GO TO 100
    
    SAT0 = S0
    
    ERR = DMDL%RUNMODEL(QFLOW, ZMAX, SAT0, TMAX, NT, TOUT, ZOUT, SATOUT, QINOUT, QOUTOUT)
    
    DEALLOCATE (SAT0)
100 RETURN
    END FUNCTION RUN_RICHARDSUC
    
    INTEGER FUNCTION RUN_RICHARDSUI(KS, L, A, N, M, DTHETA, QIN, TON, TOFF, ZMAX, NZ, TMAX, NT, S0, &
        TOUT, ZOUT, SATOUT, QINOUT, QOUTOUT) RESULT(ERR)
    !DEC$ ATTRIBUTES DLLEXPORT, ALIAS : 'RUN_RICHARDSUI' :: RUN_RICHARDSUI
    USE GENERAL
    USE INLETFLOW_MODULE
    USE UNSATURATED
    USE DARCY_MODULE
    IMPLICIT NONE
    
    REAL(DP), INTENT(IN) :: KS, L, A, N, M, DTHETA, QIN, TON, TOFF, ZMAX, TMAX, S0
    INTEGER, INTENT(IN) :: NZ, NT
    REAL(DP), INTENT(OUT) :: TOUT(NT), ZOUT(NZ), SATOUT(NZ, NT), QINOUT(NT), QOUTOUT(NT)
    
    TYPE(VANGENUCHTEN) :: UMDL
    TYPE(DARCY) :: DMDL
    TYPE(INTERMITTENT_INLETFLOW) :: QFLOW
    REAL(DP), DIMENSION(:), ALLOCATABLE :: SAT0
    
    UMDL = MKVANGENUCHTEN(L, A, N, M)
    DMDL = MKDARCY(NZ, KS, DTHETA, UMDL)
    QFLOW = INTERMITTENT_INLETFLOW(QIN, TON, TOFF)
    ERR = 0
    
    ALLOCATE (SAT0(NZ), STAT = ERR)
    IF (ERR > 0) GO TO 100
    
    SAT0 = S0
    
    ERR = DMDL%RUNMODEL(QFLOW, ZMAX, SAT0, TMAX, NT, TOUT, ZOUT, SATOUT, QINOUT, QOUTOUT)
    
    DEALLOCATE (SAT0)
100 RETURN
    END FUNCTION RUN_RICHARDSUI
    
    
    
    
    INTEGER FUNCTION RUN_RICHARDSNN(PARIN, NZ, NT, PARNOISE, TOUT, ZOUT, SATOUT, QINOUT, QOUTOUT) RESULT(ERR)
    !DEC$ ATTRIBUTES DLLEXPORT, ALIAS : 'RUN_RICHARDSNN' :: RUN_RICHARDSNN
    USE GENERAL
    USE INLETFLOW_MODULE
    USE UNSATURATED
    USE DARCY_MODULE
    IMPLICIT NONE
    
    REAL(DP), INTENT(IN), DIMENSION(11) :: PARIN
    REAL(DP), INTENT(IN), DIMENSION(6) :: PARNOISE
    INTEGER, INTENT(IN) :: NZ, NT
    REAL(DP), INTENT(OUT) :: TOUT(NT), ZOUT(NZ), SATOUT(NZ, NT), QINOUT(NT), QOUTOUT(NT)
    
    REAL(DP) :: KS, L, A, N, M, DTHETA, QIN, QNPERIOD, ZMAX, TMAX, S0
    REAL(DP) :: KSN, AN, NN, MN, QINN, S0N
    REAL(DP) :: AA, NA, MA

    TYPE(VANGENUCHTEN), DIMENSION(:), ALLOCATABLE :: UMDL
    TYPE(DARCY) :: DMDL
    TYPE(NOISY_INLETFLOW) :: QFLOW
    REAL(DP), DIMENSION(:), ALLOCATABLE :: SAT0, SKOND, DTHETA_VEC
    INTEGER :: I, NPERIODS
    
    KS = PARIN(1); L = PARIN(2); A = PARIN(3); N = PARIN(4); M = PARIN(5)
    DTHETA = PARIN(6); QIN = PARIN(7); QNPERIOD = PARIN(8); ZMAX = PARIN(9)
    TMAX = PARIN(10); S0 = PARIN(11)
    
    KSN = PARNOISE(1); AN = PARNOISE(2); NN = PARNOISE(3)
    MN = PARNOISE(4); QINN = PARNOISE(5); S0N = PARNOISE(6)
    
    ALLOCATE (UMDL(NZ+1), SAT0(NZ), SKOND(NZ+1), DTHETA_VEC(NZ), STAT=ERR)
    IF (ERR > 0) GO TO 100
    
    DTHETA_VEC = DTHETA
    
    CALL RANDOM_SEED
    
    DO I=1, NZ
        AA = NOISYVAL(A, AN)
        NA = NOISYVAL(N, NN)
        MA = NOISYVAL(M, MN)
        SAT0(I) = NOISYVAL(S0, S0N)
        ! Passing the lg of Ks, needed to make Ks
        ! vary in orders of magnitude with the noise
        SKOND(I) = 10D0 ** NOISYVAL(KS, KSN)
        UMDL(I) = MKVANGENUCHTEN(L, AA, NA, MA)
    END DO
    AA = NOISYVAL(A, AN)
    NA = NOISYVAL(N, NN)
    MA = NOISYVAL(M, MN)
    SKOND(NZ+1) = 10D0 ** NOISYVAL(KS, KSN)
    UMDL(NZ+1) = MKVANGENUCHTEN(L, AA, NA, MA)
    
    DMDL%SKONDF = SKOND
    DMDL%UMDLF = UMDL
    DMDL%DTHETAC = DTHETA_VEC
    
    NPERIODS = INT(TMAX / QNPERIOD) + 1
    QFLOW = MK_NOISY_INLETFLOW(QIN, QINN, QNPERIOD, NPERIODS)
    
    ERR = 0
    
    ERR = DMDL%RUNMODEL(QFLOW, ZMAX, SAT0, TMAX, NT, TOUT, ZOUT, SATOUT, QINOUT, QOUTOUT)
    
    DEALLOCATE (SAT0)
100 RETURN
    
    CONTAINS
    
        REAL(DP) FUNCTION NOISYVAL(BASE, NOISE)
            REAL(DP), INTENT(IN) :: BASE, NOISE
        
            REAL(DP) :: RND
        
            CALL RANDOM_NUMBER(RND)
        
            RND = 2D0*(RND-0.5D0)
        
            NOISYVAL = BASE * (1D0 + NOISE * RND)
        END FUNCTION
    
    END FUNCTION RUN_RICHARDSNN

END MODULE