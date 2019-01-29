MODULE UNSATURATED
    
    USE GENERAL
    IMPLICIT NONE
    
    TYPE, ABSTRACT, PUBLIC :: UMODEL
    CONTAINS
        PROCEDURE(IAVERAGE), PASS(THIS), DEFERRED, PUBLIC :: AVERAGE
        PROCEDURE(ICPRES), PASS(THIS), DEFERRED, PUBLIC :: CPRES
        PROCEDURE(IRCOND), PASS(THIS), DEFERRED, PUBLIC :: RCOND
        PROCEDURE(ISAT), PASS(THIS), DEFERRED, PUBLIC :: SAT
    END TYPE UMODEL
    
    
    ABSTRACT INTERFACE
        
        FUNCTION IAVERAGE(THIS, MDL2)
        IMPORT :: UMODEL
        CLASS(UMODEL), INTENT(IN) :: THIS
        CLASS(UMODEL), INTENT(IN) :: MDL2
        CLASS(UMODEL), POINTER :: IAVERAGE
        END FUNCTION IAVERAGE
        
        ELEMENTAL REAL(DP) FUNCTION ICPRES(THIS, SAT)
        IMPORT :: UMODEL, DP
        CLASS(UMODEL), INTENT(IN) :: THIS
        REAL(DP), INTENT(IN) :: SAT
        END FUNCTION ICPRES
        
        ELEMENTAL REAL(DP) FUNCTION IRCOND(THIS, SAT)
        IMPORT :: UMODEL, DP
        CLASS(UMODEL), INTENT(IN) :: THIS
        REAL(DP), INTENT(IN) :: SAT
        END FUNCTION IRCOND
        
        ELEMENTAL REAL(DP) FUNCTION ISAT(THIS, PRES)
        IMPORT :: UMODEL, DP
        CLASS(UMODEL), INTENT(IN) :: THIS
        REAL(DP), INTENT(IN) :: PRES
        END FUNCTION ISAT
    END INTERFACE
        
    
        
    TYPE, EXTENDS(UMODEL), PUBLIC :: VANGENUCHTEN
        REAL(DP) :: L, A, N, M
    CONTAINS
        PROCEDURE, PASS(THIS), PUBLIC :: AVERAGE => VG_AVERAGE
        PROCEDURE, PASS(THIS), PUBLIC :: CPRES => VG_CPRES
        PROCEDURE, PASS(THIS), PUBLIC :: RCOND => VG_RCOND
        PROCEDURE, PASS(THIS), PUBLIC :: SAT => VG_SAT
    END TYPE VANGENUCHTEN
    
    
    ! RCOND
    !   Generic subroutine, calculates relative conductivity from saturation
    !
    ! Arguments
    !   SAT         Vector of saturation (in)
    !   MDL         Given unsaturated model (in)
    !   COND        Calculated conductivity (out)
    !
    !INTERFACE RCOND
    !    MODULE PROCEDURE VG_RCOND
    !END INTERFACE RCOND
    !
    
    ! CPRES
    !   Generic subroutine, calculates capillary pressure from saturation
    !
    ! Arguments
    !   SAT         Vector of saturation (in)
    !   MDL         Given unsaturated model (in)
    !   COND        Calculated conductivity (out)
    !
    !INTERFACE CPRES
    !    MODULE PROCEDURE VG_CPRES
    !END INTERFACE CPRES
    !
    
    ! SAT
    !   Generic subroutine calculates saturation from capillary pressure
    !
    ! Arguments
    !   PRES        Vector of capillary pressure (in)
    !   MDL         Given unsaturated model (in)
    !   S           Calculated saturation (out)
    !
    !INTERFACE SAT
    !    MODULE PROCEDURE VG_SAT
    !END INTERFACE
    
    
    ! MKVANGENUCHTEN
    !   Generic van Genuchten model constructor (factory)
    !
    ! Description
    !   Convenience function to construct van Genuchten model.
    !   Has three forms:
    !   1.
    !       MDL = MKVANGENUCHTEN(L, A, N[, M])
    !   with L, A, N, M being REAL(DP) - constructs 1-dimensional model
    !   2.
    !       MDL = MKVANGENUCHTEN(I, L, A, N[, M])
    !   the same as 1, but additionally I is an INTEGER defining the
    !   dimensionality of the model.
    !   3.
    !       MDL = MKVANGENUCHTEN(VL, VA, VN, VM)
    !   with VL, VA, VN, VM being vectors of REAL(DP) - constructs
    !   van Genuchten model with given vector parameters.
    !
    !   If the form 1 or 2 is used, DELUMODEL must be used at the end of
    !   computation.
    INTERFACE MKVANGENUCHTEN
        MODULE PROCEDURE MKVANGENUCHTEN1, MKVANGENUCHTENN, MKVANGENUCHTENV
    END INTERFACE
    
    
    ! DELMODEL
    !   Deallocates unsaturated model parameters. Whether it must be used on
    !   a particular model, see the documentation of the model constructor.
    INTERFACE DELUMODEL
        MODULE PROCEDURE DELVANGENUCHTEN
    END INTERFACE
    
    !INTERFACE AVERAGEUMODEL
    !    MODULE PROCEDURE AVERAGEUMODEL_VG
    !END INTERFACE
    
    CONTAINS
    
    
    ELEMENTAL REAL(DP) FUNCTION VG_RCOND(THIS, SAT)
    ! VG_RCOND
    !   Calculates relative conductivity for the van Genuchten model
    !
    ! Arguments
    !   SAT         Input saturation
    !   MDL         (Input) model
    !   KOND        Output: relative conductivities
    !
    CLASS(VANGENUCHTEN), INTENT(IN) :: THIS
    REAL(DP), INTENT(IN) :: SAT
    
    REAL(DP) :: S
    
    IF (SAT >= 1D0) THEN
        S = 0.999D0
    ELSE IF (SAT <= 0D0) THEN
        S = 0.001D0
    ELSE
        S = SAT
    END IF

    
    VG_RCOND = (S ** THIS%L) * &
        (1.0D0 - (1.0D0 - S ** (1.0D0 / THIS%M)) ** THIS%M) ** 2
    
    END FUNCTION VG_RCOND
    
    
    ELEMENTAL REAL(DP) FUNCTION VG_CPRES(THIS, SAT)
    ! VG_CPRES
    !   Calculates capillary pressure for van Genuchten model
    !
    ! Arguments
    !   SAT         Input saturation
    !   MDL         (Input) model
    !   PRES        Output: capillary pressure
    !
    REAL(DP), INTENT(IN) :: SAT
    CLASS(VANGENUCHTEN), INTENT(IN) :: THIS
    
    REAL(DP) :: S
    
    IF (SAT >= 1D0) THEN
        S = 0.999D0
    ELSE IF (SAT <= 0D0) THEN
        S = 0.001D0
    ELSE
        S = SAT
    END IF

    
    VG_CPRES = - 1.0D0 / THIS%A * &
        (S ** (-1.0D0 / THIS%M) - 1.0D0) ** (1.0D0 / THIS%N)
    
    END FUNCTION VG_CPRES
    
    
    ELEMENTAL REAL(DP) FUNCTION VG_SAT(THIS, PRES)
    ! V_SAT
    !   Calculates saturation from pressure of van Genuchten model
    !
    ! Arguments
    !   PRES        Input capillary pressure
    !   MDL         (Input) van Genuchten model
    !   PRES        Output: capillary pressure
    !
    REAL(DP), INTENT(IN) :: PRES
    CLASS(VANGENUCHTEN), INTENT(IN) :: THIS
    
    VG_SAT = (1.0D0 / (1.0D0 + ABS(THIS%A * PRES) ** THIS%N)) ** THIS%M
    
    END FUNCTION VG_SAT
    
    
    FUNCTION MKVANGENUCHTEN1(L, A, N, M) RESULT(VG)
    ! MKVANGENUCHTEN1
    !   Constructs one-dimensional van Genuchten model.
    !
    ! Arguments
    !   L, A, N, M  Input REAL(DP) model parameters
    !                   M is optional, if not provided: M = 1-1/N
    !
    ! Result
    !   van Genuchten model (TYPE(VANGENUCHTEN))
    !
    REAL(DP), INTENT(IN) :: L, A, N
    REAL(DP), INTENT(IN), OPTIONAL :: M
    
    TYPE(VANGENUCHTEN) :: VG
    
    REAL(DP) :: M1
        
    IF (.NOT. PRESENT(M)) THEN
        M1 = 1.0D0 - 1.0D0 / N
    ELSE
        M1 = M
    END IF
    
    VG%L = L; VG%A = A; VG%N = N; VG%M = M1
    
    END FUNCTION MKVANGENUCHTEN1
    
    
    FUNCTION MKVANGENUCHTENN(I, L, A, N, M) RESULT(VG)
    ! MKVANGENUCHTENN
    !   Constructs I-dimensional van Genuchten model.
    !   It allocates space for the model in a vector - DELUMODEL
    !   must be used at the end of computation.
    !
    ! Arguments
    !   I           Model dimensionality
    !   L, A, N, M  Input REAL(DP) model parameters
    !                   M is optional, if not provided: M = 1-1/N
    !
    ! Result
    !   Vector of van Genuchten models (TYPE(VANGENUCHTEN), DIMENSION(I))
    !
    INTEGER, INTENT(IN) :: I
    REAL(DP), INTENT(IN) :: L, A, N
    REAL(DP), INTENT(IN), OPTIONAL :: M
    
    TYPE(VANGENUCHTEN), DIMENSION(:), ALLOCATABLE :: VG
    
    INTEGER :: ERR
    REAL(DP) :: M1
    
    ALLOCATE (VG(I), STAT=ERR)
    
    IF (ERR > 0) THEN
        WRITE (*,*) "MKVANGENUCHTENN: Cannot allocate space. Stopping"
        STOP
    END IF
    
    IF (.NOT. PRESENT(M)) THEN
        M1 = 1.0D0 - 1.0D0 / N
    ELSE
        M1 = M
    END IF
    
    VG(:)%L = L; VG(:)%A = A; VG(:)%N = N; VG(:)%M = M1
    
    END FUNCTION MKVANGENUCHTENN
    
    
    FUNCTION MKVANGENUCHTENV(L, A, N, M) RESULT(VG)
    ! MKVANGENUCHTENV
    !   Constructs general multi-dimensional van Genuchten model.
    !   Allocates the space for the model - use DELUMODEL at the end of
    !   computation.
    !
    ! Arguments
    !   L, A, N, M  Input REAL(DP) model parameters
    !               NOTE: M is compulsary in this form! 
    !
    ! Result
    !   Vector of van Genuchten models (TYPE(VANGENUCHTEN))
    !
    REAL(DP), DIMENSION(:), INTENT(IN) :: L, A, N
    REAL(DP), DIMENSION(:), OPTIONAL, INTENT(IN) :: M
    
    TYPE(VANGENUCHTEN), DIMENSION(:), ALLOCATABLE :: VG
    
    INTEGER :: I, ERR
    
    I = SIZE(L)
    ALLOCATE (VG(I), STAT=ERR)
    IF (ERR > 0) THEN
        WRITE (*,*) "VANGENUCHTENV: Cannot allocate space for the model. Stopping"
        STOP
    END IF
    
    VG(:)%L = L(:); VG(:)%A = A(:); VG(:)%N = N(:)
    IF (PRESENT(M)) THEN
        VG(:)%M = M(:)
    ELSE
        VG(:)%M = 1D0 - 1D0/N(:)
    END IF
    
    END FUNCTION MKVANGENUCHTENV
    
    
    SUBROUTINE DELVANGENUCHTEN(VG)
    ! DELVANGENUCHTEN
    !   Deallocates the storage of van Genuchten model
    !   The model and any of its storage is unusable after
    !   this call. Trying to access them will cause the program
    !   to crush
    !
    ! Arguments
    !   VG          Model to be deallocated
    !
    TYPE(VANGENUCHTEN), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: VG
    
    INTEGER :: ERR
    
    ERR = 0    
    ! Only deallocate if the array is allocated
    ! If error occurs - go to the error trap
    IF (ALLOCATED(VG)) DEALLOCATE (VG, STAT=ERR)
    IF (ERR > 0) GO TO 100
        
    RETURN
    ! Error trap
100 WRITE (*,*) "DELVANGENUCHTEN: Cannot deallocate the model. Stopping"
    STOP
    
    END SUBROUTINE DELVANGENUCHTEN
    
    
    !FUNCTION UMODEL_AVERAGE(MDL1, MDL2) RESULT(MDL)
    !CLASS(UMODEL), INTENT(IN) :: MDL1, MDL2
    !TYPE(UMODEL), TARGET :: M
    !CLASS(UMODEL), POINTER :: MDL
    !MDL => M
    !END FUNCTION
    
    FUNCTION VG_AVERAGE(THIS, MDL2) RESULT(MDL)
    ! AVERAGEVG
    !   Constructs the arithmetic average of two models
    !   Purpose: to construct centroid models based on face-models
    !           in finite-volume discretization
    CLASS(VANGENUCHTEN), INTENT(IN) :: THIS
    CLASS(UMODEL), INTENT(IN) :: MDL2
    CLASS(UMODEL), POINTER :: MDL
    TYPE(VANGENUCHTEN) :: MDL2V
    TYPE(VANGENUCHTEN), TARGET :: MDLV
    
    MDL2V = TRANSFER(MDL2, MOLD=THIS)
    MDLV%L = 0.5D0*(THIS%L + MDL2V%L)
    MDLV%A = 0.5D0*(THIS%A + MDL2V%A)
    MDLV%N = 0.5D0*(THIS%N + MDL2V%N)
    MDLV%M = 0.5D0*(THIS%M + MDL2V%M)
    
    MDL => MDLV
    
    END FUNCTION VG_AVERAGE
    
    !ELEMENTAL FUNCTION AVERAGEUMODEL(MDL1, MDL2) RESULT(MDL)
    !! AVERAGEUMODEL
    !!   Averages to unsaturated models by means of arithmetic average
    !!   Example of use:
    !!       AVUMDL = AVERAGEUMODEL(MDL1, MDL2)
    !!   This function is required to be elemental.
    !CLASS(UMODEL), INTENT(IN) :: MDL1, MDL2
    !CLASS(UMODEL), POINTER :: MDL
    !MDL = MDL1 % AVERAGE(MDL2)
    !END FUNCTION
    
    
END MODULE UNSATURATED