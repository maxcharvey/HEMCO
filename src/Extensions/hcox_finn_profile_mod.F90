! Opt-in GFAS shape transfer; FINNv2.5 remains the sole mass inventory.
MODULE HCOX_FINN_PROFILE_MOD
  USE HCO_ERROR_MOD
  USE HCO_STATE_MOD, ONLY: HCO_State
  USE HCOX_STATE_MOD, ONLY: Ext_State
  USE HCOX_FINN_PROFILE_KERNEL_MOD, ONLY: FINN_Normalize, FINN_Allocate
  USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_FINITE
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: FINN_Profile_Init, FINN_Profile_Prepare, FINN_Profile_Apply
  PUBLIC :: FINN_Profile_Report, FINN_Profile_Final, ProfileEnabled
  LOGICAL, SAVE :: ProfileEnabled = .FALSE.
  CHARACTER(LEN=31), PARAMETER :: MethodName = 'gfas_prepared'
  CHARACTER(LEN=31), SAVE :: Fallback = 'pbl'
  REAL(hp), ALLOCATABLE, SAVE :: Weights(:,:,:)
  INTEGER, ALLOCATABLE, SAVE :: Reason(:,:) ! 0 valid, 1 unavailable, 2 terrain
CONTAINS
  SUBROUTINE FINN_Profile_Init( HcoState, ExtState, ExtNr, IDs, Names, RC )
    USE HCO_EXTLIST_MOD, ONLY: GetExtOpt
    USE HCO_DIAGN_MOD, ONLY: Diagn_Create
    TYPE(HCO_State), POINTER :: HcoState
    TYPE(Ext_State), POINTER :: ExtState
    INTEGER, INTENT(IN) :: ExtNr, IDs(:)
    CHARACTER(LEN=*), INTENT(IN) :: Names(:)
    INTEGER, INTENT(INOUT) :: RC
    LOGICAL :: Found
    INTEGER :: N, M, AS
    CHARACTER(LEN=20), PARAMETER :: Metrics(4) = [ CHARACTER(LEN=20) :: &
         'FINNProfileSource_', 'FINNProfileFallback_', &
         'FINNProfileTerrain_', 'FINNProfileAbovePBL_' ]
    CALL GetExtOpt(HcoState%Config, ExtNr, 'FINNV25_GFAS_PROFILE', &
                   OptValBool=ProfileEnabled, FOUND=Found, RC=RC)
    IF ( RC /= HCO_SUCCESS ) RETURN
    IF ( .NOT. Found ) ProfileEnabled = .FALSE.
    IF ( .NOT. ProfileEnabled ) RETURN
    CALL GetExtOpt(HcoState%Config, ExtNr, 'FINNv25_profile_fallback', &
                   OptValChar=Fallback, FOUND=Found, RC=RC)
    IF ( RC /= HCO_SUCCESS ) RETURN
    IF ( .NOT. Found ) Fallback = 'pbl'
    IF ( Fallback /= 'pbl' .AND. Fallback /= 'surface' .AND. Fallback /= 'error' ) THEN
       CALL HCO_ERROR('FINN profile fallback must be pbl, surface, or error', RC)
       RETURN
    ENDIF
    ALLOCATE(Weights(HcoState%NX,HcoState%NY,HcoState%NZ), &
             Reason(HcoState%NX,HcoState%NY), STAT=AS)
    IF ( AS /= 0 ) THEN
       CALL HCO_ERROR('Cannot allocate FINN profile arrays', RC)
       RETURN
    ENDIF
    ! Needed also for the emitted above-PBL metric, in every hybrid mode.
    ExtState%FRAC_OF_PBL%DoUse = .TRUE.
    ExtState%PBL_OCCUPANCY%DoUse = .TRUE.
    DO N=1,SIZE(IDs)
       IF ( IDs(N) < 0 ) CYCLE
       DO M=1,SIZE(Metrics)
          CALL Diagn_Create(HcoState=HcoState, cName=TRIM(Metrics(M))//TRIM(Names(N)), &
               ExtNr=ExtNr, Cat=-1, Hier=-1, HcoID=IDs(N), SpaceDim=2, &
               OutUnit='kg/m2/s', AutoFill=0, RC=RC)
          IF ( RC /= HCO_SUCCESS ) RETURN
       ENDDO
    ENDDO
    IF ( HcoState%amIRoot ) THEN
       CALL HCO_MSG('FINNv2.5 auxiliary profile: '//TRIM(MethodName), LUN=HcoState%Config%hcoLogLUN)
       CALL HCO_MSG('FINNv2.5 profile fallback: '//TRIM(Fallback), LUN=HcoState%Config%hcoLogLUN)
    ENDIF
  END SUBROUTINE FINN_Profile_Init

  SUBROUTINE FINN_Profile_Prepare( HcoState, RC )
    USE HCO_CALC_MOD, ONLY: HCO_EvalFld
    TYPE(HCO_State), POINTER :: HcoState
    INTEGER, INTENT(INOUT) :: RC
    REAL(hp), ALLOCATABLE :: Ref(:,:,:)
    INTEGER :: I, J, S, AS
    ALLOCATE(Ref(HcoState%NX,HcoState%NY,HcoState%NZ), STAT=AS)
    IF ( AS /= 0 ) THEN
       CALL HCO_ERROR('Cannot allocate FINN auxiliary reference', RC)
       RETURN
    ENDIF
    Ref=0.0_hp
    CALL HCO_EvalFld(HcoState, 'FINNV25_GFAS_REFERENCE', Ref, RC)
    IF ( RC /= HCO_SUCCESS ) RETURN
    Reason=0
    DO J=1,HcoState%NY
    DO I=1,HcoState%NX
       CALL FINN_Normalize(Ref(I,J,:), Weights(I,J,:), S)
       IF ( S == 1 ) THEN
          Reason(I,J)=1
       ELSEIF ( S /= 0 ) THEN
          CALL HCO_ERROR('Invalid FINN GFAS reference: require finite nonnegative values', RC)
          RETURN
       ENDIF
    ENDDO
    ENDDO
  END SUBROUTINE FINN_Profile_Prepare

  SUBROUTINE FINN_Profile_Apply( HcoState, ExtState, Source, Emission, RC )
    TYPE(HCO_State), POINTER :: HcoState
    TYPE(Ext_State), POINTER :: ExtState
    REAL(hp), INTENT(IN) :: Source(:,:)
    REAL(hp), INTENT(OUT) :: Emission(:,:,:)
    INTEGER, INTENT(INOUT) :: RC
    INTEGER :: I, J, S
    REAL(hp) :: W(HcoState%NZ), PBL(HcoState%NZ)
    Emission=0.0_hp
    IF ( ANY(.NOT. IEEE_IS_FINITE(Source)) .OR. ANY(Source < 0.0_hp) ) THEN
       CALL HCO_ERROR('Invalid FINNv2.5 source for profile allocation', RC)
       RETURN
    ENDIF
    DO J=1,HcoState%NY
    DO I=1,HcoState%NX
       IF ( Source(I,J) == 0.0_hp ) CYCLE
       PBL=ExtState%FRAC_OF_PBL%Arr%Val(I,J,:)
       IF ( ANY(.NOT. IEEE_IS_FINITE(PBL)) .OR. ANY(PBL<0.0_hp) .OR. ANY(PBL>1.0_hp) ) THEN
          CALL HCO_ERROR('Invalid PBL fractions in FINN profile allocation', RC)
          RETURN
       ENDIF
       PBL=ExtState%PBL_OCCUPANCY%Arr%Val(I,J,:)
       IF ( ANY(.NOT. IEEE_IS_FINITE(PBL)) .OR. ANY(PBL<0.0_hp) .OR. ANY(PBL>1.0_hp) ) THEN
          CALL HCO_ERROR('Invalid native PBL occupancy in FINN profile diagnostics', RC)
          RETURN
       ENDIF
       PBL=ExtState%FRAC_OF_PBL%Arr%Val(I,J,:)
       W=Weights(I,J,:)
       IF ( Reason(I,J) /= 0 ) THEN
          SELECT CASE (TRIM(Fallback))
          CASE ('pbl')
             CALL FINN_Normalize(PBL,W,S)
             IF ( S /= 0 ) THEN
                CALL HCO_ERROR('FINN profile PBL fallback has no valid PBL support', RC)
                RETURN
             ENDIF
          CASE ('surface')
             W=0.0_hp
             W(1)=1.0_hp
          CASE DEFAULT
             CALL HCO_ERROR('FINN emissions lack GFAS profile support; fallback=error', RC)
             RETURN
          END SELECT
       ENDIF
       CALL FINN_Allocate(Source(I,J),W,Emission(I,J,:),S)
       IF ( S /= 0 ) THEN
          CALL HCO_ERROR('FINN profile allocation failed positivity or closure', RC)
          RETURN
       ENDIF
    ENDDO
    ENDDO
  END SUBROUTINE FINN_Profile_Apply

  SUBROUTINE FINN_Profile_Report( HcoState, ExtState, ID, Name, Source, Emission, RC )
    USE HCO_DIAGN_MOD, ONLY: Diagn_Update
    USE HCO_SCALE_MOD, ONLY: HCO_ScaleArr
    TYPE(HCO_State), POINTER :: HcoState
    TYPE(Ext_State), POINTER :: ExtState
    INTEGER, INTENT(IN) :: ID
    CHARACTER(LEN=*), INTENT(IN) :: Name
    REAL(hp), INTENT(IN) :: Source(:,:), Emission(:,:,:)
    INTEGER, INTENT(INOUT) :: RC
    REAL(hp) :: A(HcoState%NX,HcoState%NY), Column(HcoState%NX,HcoState%NY)
    INTEGER :: I,J
    ! Called after EmisAdd: Emission already includes universal species scaling.
    ! Independently scale a COPY of input for the source closure diagnostic.
    A=Source
    CALL HCO_ScaleArr(HcoState,ID,A,RC)
    IF ( RC /= HCO_SUCCESS ) RETURN
    CALL Diagn_Update(HcoState,cName='FINNProfileSource_'//TRIM(Name), &
                       AutoFill=0,Array2D=A,RC=RC)
    IF ( RC /= HCO_SUCCESS ) RETURN
    Column=SUM(Emission,DIM=3)
    A=0.0_hp
    WHERE (Reason/=0) A=Column
    CALL Diagn_Update(HcoState,cName='FINNProfileFallback_'//TRIM(Name), &
                       AutoFill=0,Array2D=A,RC=RC)
    IF ( RC /= HCO_SUCCESS ) RETURN
    A=0.0_hp
    WHERE (Reason==2) A=Column
    CALL Diagn_Update(HcoState,cName='FINNProfileTerrain_'//TRIM(Name), &
                       AutoFill=0,Array2D=A,RC=RC)
    IF ( RC /= HCO_SUCCESS ) RETURN
    A=0.0_hp
    DO J=1,HcoState%NY
    DO I=1,HcoState%NX
       IF ( Source(I,J)>0.0_hp ) &
          A(I,J)=SUM(Emission(I,J,:)*(1.0_hp-ExtState%PBL_OCCUPANCY%Arr%Val(I,J,:)))
    ENDDO
    ENDDO
    CALL Diagn_Update(HcoState,cName='FINNProfileAbovePBL_'//TRIM(Name), &
                       AutoFill=0,Array2D=A,RC=RC)
  END SUBROUTINE FINN_Profile_Report

  SUBROUTINE FINN_Profile_Final()
    IF ( ALLOCATED(Weights) ) DEALLOCATE(Weights)
    IF ( ALLOCATED(Reason) ) DEALLOCATE(Reason)
    ProfileEnabled=.FALSE.
    Fallback='pbl'
  END SUBROUTINE FINN_Profile_Final
END MODULE HCOX_FINN_PROFILE_MOD
