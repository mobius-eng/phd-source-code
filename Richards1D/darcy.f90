MODULE DARCY_MODULE
    
    USE GENERAL
    USE UNSATURATED
    USE INLETFLOW_MODULE
    IMPLICIT NONE
    
    
    ! DARCY
    TYPE DARCY
        REAL(DP), DIMENSION(:), ALLOCATABLE :: SKONDF, DTHETAC
        CLASS(UMODEL), DIMENSION(:), ALLOCATABLE :: UMDLF
    CONTAINS
        PROCEDURE :: RUNMODEL => DARCY_RUNMODEL
    END TYPE DARCY
    
    
    INTERFACE MKDARCY
        MODULE PROCEDURE MKDARCY_UNIFORM
    END INTERFACE
    
    
    CONTAINS
    
    
    TYPE(DARCY) FUNCTION MKDARCY_UNIFORM(NZC, SKOND, DTHETA, UMDL) RESULT(DMDL)
    CLASS(UMODEL), INTENT(IN) :: UMDL
    REAL(DP), INTENT(IN) :: SKOND, DTHETA
    INTEGER :: NZC, ERR
    
    ALLOCATE(DMDL%SKONDF(NZC+1), DMDL%DTHETAC(NZC),STAT=ERR)
    IF (ERR > 0) GO TO 100
    ALLOCATE(DMDL%UMDLF(NZC+1), MOLD=UMDL, STAT=ERR)
    IF (ERR > 0) GO TO 100
    
    DMDL%SKONDF = SKOND
    DMDL%DTHETAC = DTHETA
    DMDL%UMDLF = UMDL
    RETURN
    
100 WRITE (*,*) "ERROR: Cannot allocate the space for the Darcy model. Stopping"
    STOP
    END FUNCTION MKDARCY_UNIFORM
    

    
    
    SUBROUTINE DELDARCY(DMDL)
    TYPE(DARCY), INTENT(INOUT) :: DMDL
    
    INTEGER :: ERR
    
    DEALLOCATE (DMDL%SKONDF, DMDL%DTHETAC, DMDL%UMDLF, STAT=ERR)
    IF (ERR > 0) GO TO 100
    
    RETURN
    
100 WRITE (*,*) "ERROR: Could not deallocate Darcy's model"
    END SUBROUTINE DELDARCY
    
    
    
    INTEGER FUNCTION DARCY_RUNMODEL(DMDL, QIN, ZMAX, SAT0, TF, NT, &
        TOUT, ZOUT, SATOUT, QINOUT, QOUTOUT) RESULT(IRESULT)
    ! DARCY_RUNMODEL
    !   The driver of the model simulation: calculates the evolution of
    !   saturation with space and time.
    !   It is a type-bound procedure of type DARCY.
    !
    ! Arguments
    !       DMDL            (IN) Type DARCY instance
    !       QIN             (IN) Type INLETFLOW (or derivative) instance
    !       ZMAX            (IN) Maximum spatial depth (height)
    !       SAT0            (IN) Vector of initial saturation
    !       TF              (IN) Final time of simulation
    !       NT              (IN) Number of output time points requrired
    !       TOUT            (OUT) Vector of time points of simulation
    !       ZOUT            (OUT) Vector of spatial points (centroids)
    !       SATOUT          (OUT) Matrix of calculated saturation
    !       QINOUT          (OUT) Vector of stored inlet flow
    !       QOUTOUT         (OUT) Vector of stored outlet flow
    !
    ! Notes
    !   The sizes of all arrays must be conforming:
    !       NT = LENGTH(TOUT), NZC = LENGTH(SAT0) = LENGTH(ZOUT)
    !       SIZE(SATOUT) = (NT, NZC)
    !       NT = LENGTH(QINOUT) = LENGTH(QOUTOUT)
    !   Also, the sizes of arrays in DMDL must be:
    !       DMDL%UMDLF, DMDL%SKONDF     NZC+1
    !       DMDL%DTHETAC                NZC
    ! -------------     Block 0: variable declarations      -------------------
    INTERFACE
        SUBROUTINE DLSODA(F, NEQ, Y, T, TOUT, ITOL, RTOL, ATOL, ITASK, &
            ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JAC, JT)
        USE GENERAL
        EXTERNAL :: F, JAC
        INTEGER :: NEQ, ITOL, ITASK, ISTATE, IOPT, LRW, LIW, JT, IWORK
        REAL(DP) :: T, TOUT, RTOL, Y, ATOL, RWORK
        DIMENSION Y(*), RWORK(LRW), IWORK(LIW)
        END SUBROUTINE DLSODA
    END INTERFACE
    ! Input arguments
    CLASS(DARCY), INTENT(IN) :: DMDL                ! DARCY object
    CLASS(INLETFLOW), INTENT(INOUT) :: QIN             ! Inlet flow object
    REAL(DP), INTENT(IN) :: ZMAX, TF                ! Max depth and time
    REAL(DP), DIMENSION(:), INTENT(IN) :: SAT0      ! Initial saturation
    INTEGER, INTENT(IN) :: NT                       ! Number of time outputs required
    ! Outputs
    REAL(DP), DIMENSION(:), INTENT(OUT) :: TOUT, &  ! Arrays of time and space points,
        ZOUT, QINOUT, QOUTOUT                       ! in and out flow with time
    REAL(DP), DIMENSION(:,:), INTENT(OUT) :: SATOUT ! matrix of saturaion (time, space)
    ! Internal variables
    !   NZC     Number of centroids
    !   NZF(I)  Number of (internal) faces
    !   ERR     Error indicator
    !   ITCUR   Time-step counter
    !   I       General counter
    !   DZ      Size (height) of the volume
    !   JF      Flux on faces
    !   SATAVFI Average saturation on internal faces
    !   KONDFI  Conductivity on internal faces
    !   PRESC   Capillary pressure on centroids
    !   UMDLC   Unsaturated models on centroids (average between faces)
    INTEGER :: NZC, NZF, NZFI, ERR, ITCUR, I, ITRYS
    REAL(DP) :: DZ
    REAL(DP), DIMENSION(:), ALLOCATABLE :: JF, SATAVFI, KONDFI, PRESC
    CLASS(UMODEL), DIMENSION(:), ALLOCATABLE :: UMDLC
    ! LSODA parameters
    INTEGER, DIMENSION(:), ALLOCATABLE ::  IWORK
    INTEGER :: NEQ, ITOL, ITASK, ISTATE, IOPT, LRW, LIW, JT
    REAL(DP) :: TCUR, TNEXT, RTOL, ATOL
    REAL(DP), DIMENSION(:), ALLOCATABLE :: SATCCUR, RWORK
    ! Dummy Jacobian function
    EXTERNAL JDUM
    ! -----------               End of Block 0             --------------------
    ! -----------            Block A: sanity checks        --------------------
    IF (.NOT. (NT == SIZE(TOUT) .AND. NT == SIZE(QINOUT) &
        .AND. NT == SIZE(QOUTOUT) .AND. NT == SIZE(SATOUT, 2))) &
        GO TO 200
    
    IF (.NOT. (SIZE(SAT0) == SIZE(SATOUT, 1) .AND. &
        SIZE(SAT0) == SIZE(DMDL%DTHETAC))) GO TO 210
    
    IF (.NOT. (ALL(SAT0 > 0D0) .AND. ALL(SAT0 < 1D0))) GO TO 220
    
    IF (TF <= 0D0) GO TO 230
    
    IF (ZMAX <= 0D0) GO TO 240
    
    IF (NT <= 0) GO TO 250
    ! -----------               End of Block A             --------------------
    ! -----------       Block B: initialize parameters     --------------------
    NZC = SIZE(DMDL%DTHETAC)                ! Number of centroids
    NZF = NZC + 1                           ! Number of faces
    NZFI = NZC - 1                          ! Number of internal faces
    DZ = ZMAX / NZC                         ! One element height
    ! Allocate internal storage
    ALLOCATE (JF(NZF), SATAVFI(NZFI), KONDFI(NZFI), PRESC(NZC), &
        STAT=ERR)
    ALLOCATE(UMDLC(NZC),MOLD=DMDL%UMDLF(1))
    IF (ERR > 0) GO TO 100
    ! Fill in TOUT and ZOUT
    CALL LINSPACE(0D0, TF, NT, TOUT)
    CALL LINSPACE(0.5D0*DZ, ZMAX-0.5D0*DZ, NZC, ZOUT)
    ! Fill unsaturated model on centroids
    FORALL (I=1:NZC) UMDLC(I) = DMDL%UMDLF(I)%AVERAGE(DMDL%UMDLF(I+1))
    ! Allocate storage for LSODA
    ALLOCATE (IWORK(20+NZC), SATCCUR(NZC), &
        RWORK(100 + 22 + NZC * MAX(16, NZC + 9)), STAT=ERR)
    IF (ERR > 0) GO TO 100
    ITRYS = 0
    ! Init LSODA data
    RWORK(:) = -5D0                         ! Useful for debugging
    IWORK(:) = 15                           ! Useful for debugging
    IWORK(1) = 1                            ! Upper and lower bands
    IWORK(2) = 1                            !   in Jacobian
    NEQ = NZC                               ! Number of equations
    RTOL = 1.0D-5                           ! Relative tolerance
    ATOL = 1.0D-5                           ! Absolute tolerance
    ITOL = 1                                ! ATOL is scalar
    ITASK = 1                               ! Initial task, always 1
    ISTATE = 1                              ! Initial state
    IOPT = 0                                ! No options
    LRW = 22 + 100 + NZC * MAX(16, NZC + 9) ! Size of RWORK
    LIW = 20 + NZC                          ! Size of IWORK
    JT = 5                                  ! No Jacobian, use banded
    SATCCUR(:) = SAT0(:)                    ! Set initial value
    SATOUT(:,1) = SAT0(:)                   ! Store initial value
    QINOUT(1) = QIN % Q(0D0)                ! The first value of flow in
    ! ---------------          End of Block B         -------------------------
    ! ---------------       Block C: solution of ODE  -------------------------
    DO ITCUR=1,NT-1
        TCUR = TOUT(ITCUR)
        TNEXT = TOUT(ITCUR+1)
10      CALL DLSODA(FDARCY, NEQ, SATCCUR, TCUR, TNEXT, ITOL, RTOL, ATOL, ITASK, &
            ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JDUM, JT)
        IF (ISTATE == -1 .AND. ITRYS < 20) THEN
            WRITE (*,20) "DLSODA: too much work encountered: Restarting", TCUR, TNEXT
20          FORMAT(/, A, 'TCUR = ', G8.4, 4X,'TCUR = ',G8.4) 
            ISTATE = 2
            ITRYS = ITRYS + 1
            GO TO 10
        END IF
        IF (ISTATE < 0) GO TO 110           ! Error has occurred
        SATOUT(:, ITCUR+1) = SATCCUR(:)     ! Record saturation
        QINOUT(ITCUR+1) = QIN % Q(TNEXT)    ! Record inlet flow
        WRITE (*,'(A)', ADVANCE='NO') "="   ! Make a mark
        !20      FORMAT(' At T =',D8.0,'   S(2) =',D14.6)
    END DO
    ! Set the outlet flow = conductivity at the bottom
    QOUTOUT(:) = DMDL%UMDLF(NZF)%RCOND(SATOUT(NZC,:))
    QOUTOUT(:) = DMDL%SKONDF(NZF) * QOUTOUT(:)
    ! ---------------------       End of Block C     --------------------------
    ! ---------------------   Block D: clearing up   --------------------------
    DEALLOCATE (JF, SATAVFI, KONDFI, PRESC, UMDLC, IWORK, RWORK, STAT=ERR)
    IF (ERR > 0) GO TO 110
    
    IRESULT = 0
    RETURN
    ! -------------            End of Block D         -------------------------
    ! -------------         Block E: Error traps      -------------------------
100 WRITE (*,*) "ERROR in DARCY_RUNMODEL: cannot allocate work space"
    IRESULT = -1 - 10
    RETURN
    
110 WRITE (*,*) "ERROR in DARCY_RUNMODEL: error in LSODA", ISTATE
    IRESULT = ISTATE - 100
    RETURN
    
120 WRITE (*,*) "ERROR in DARCY_RUNMODEL: cannot deallocate the work space"
    IRESULT = -2 - 10
    RETURN
    
200 WRITE (*,*) "ERROR in DARCY_RUNMODEL: time-related array sizes do not conform"
    IRESULT = -1 - 200
    RETURN
    
210 WRITE (*,*) "ERROR in DARCY_RUNMODEL: space-related array sizes do not conform"
    IRESULT = -2 - 200
    RETURN
    
220 WRITE (*,*) "ERROR in DARCY_RUNMODEL: initial saturation values are not in (0,1)"
    IRESULT = -3 - 200
    RETURN
    
230 WRITE (*,*) "ERROR in DARCY_RUNMODEL: final time TF is non-positive"
    IRESULT = -4 - 200
    
240 WRITE (*,*) "ERROR in DARCY_RUNMODEL: depth ZMAX is not positive"
    IRESULT = -5 - 200
    RETURN
    
250 WRITE (*,*) "ERROR in DARCY_RUNMODEL: number of time steps is not positive"
    IRESULT = -6 - 200
    RETURN
    ! -------------           End of Block E            -----------------------
    
    ! -------------   Block F: internal procedures      -----------------------
    CONTAINS
    
    SUBROUTINE FDARCY(NEQ, T, SATC, SATCDOT)
    ! FDARCY
    !   A callback function for ODE solver: RHS of ODE
    ! -------------    Block 0: variables declaration   -----------------------
    INTEGER, INTENT(IN) :: NEQ                          ! Number of equations = NZC
    REAL(DP) :: T                                       ! Current time
    REAL(DP), DIMENSION(NEQ), INTENT(IN) :: SATC        ! Current saturation
    REAL(DP), DIMENSION(NEQ), INTENT(OUT) :: SATCDOT    ! Derivative ds/dt
    INTEGER :: I                                        ! General counter
    ! ----------------         End of Block 0           -----------------------
    ! ----------------    Block A: main calculation     -----------------------
    PRESC(:) = UMDLC(:)%CPRES(SATC(:))                  ! Calculate pressure on centroids
    SATAVFI(:) = 0.50D0*(SATC(1:NZC-1) + SATC(2:NZC))   ! Av. saturation on internal faces
    KONDFI(:) = DMDL%UMDLF(2:NZF-1)%RCOND(SATAVFI)      ! Relative conductivities
    KONDFI = DMDL%SKONDF(2:NZF-1) * KONDFI              ! Full cond. on internal faces
    ! Flux
    JF(1) = QIN % Q(T)                                  ! Top volume: inlet flux
    JF(2:NZF-1) = - KONDFI * &                          ! Flux on internal faces:
        ((PRESC(2:NZC) - PRESC(1:NZC-1))/DZ - 1.0D0)    !   use Darcy's law
    JF(NZF) = DMDL%UMDLF(NZF)%RCOND(SATC(NZC))
    JF(NZF) = DMDL%SKONDF(NZF) * JF(NZF)                !   free flow
    
    SATCDOT(1:NEQ) = -1D0 / DMDL%DTHETAC(1:NZC) * (JF(2:NZF) - JF(1:NZF-1)) / DZ
    ! -----------------         End of Block A          -----------------------
    END SUBROUTINE FDARCY
    
    END FUNCTION DARCY_RUNMODEL
    
    
END MODULE DARCY_MODULE
    
SUBROUTINE JDUM()
END SUBROUTINE
    