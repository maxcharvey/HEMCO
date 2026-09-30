! Opt-in GFAS shape transfer; FINNv2.5 remains the sole mass inventory.
MODULE HCOX_FINN_PROFILE_MOD
  USE HCO_ERROR_MOD
  USE HCO_STATE_MOD, ONLY: HCO_State
  USE HCOX_STATE_MOD, ONLY: Ext_State
  USE HCOX_FINN_PROFILE_KERNEL_MOD, ONLY: FINN_Normalize, FINN_Allocate, FINN_HeightProfile
  USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_FINITE
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: FINN_Profile_Init, FINN_Profile_Prepare, FINN_Profile_Apply
  PUBLIC :: FINN_Profile_Report, FINN_Profile_Final, ProfileEnabled
  LOGICAL, SAVE :: ProfileEnabled = .FALSE.
  CHARACTER(LEN=31), PARAMETER :: MethodName = 'gfas_uniform'
  CHARACTER(LEN=31), SAVE :: Fallback = 'pbl'
  CHARACTER(LEN=31), SAVE :: HeightDatum='UNSET'
  REAL(hp), ALLOCATABLE, SAVE :: SupportFraction(:,:)
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
    CHARACTER(LEN=32), PARAMETER :: Metrics(5) = [ CHARACTER(LEN=32) :: &
         'FINNProfileSource_', 'FINNProfileFallback_', &
         'FINNProfileTerrain_', 'FINNProfileAbovePBL_', 'FINNProfileMissingSupport_' ]
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
    CALL GetExtOpt(HcoState%Config, ExtNr, 'FINNv25_profile_datum', &
                   OptValChar=HeightDatum, FOUND=Found, RC=RC)
    IF (RC/=HCO_SUCCESS) RETURN
    IF (.NOT.Found .OR. (HeightDatum/='MSL' .AND. HeightDatum/='AGL')) THEN
       CALL HCO_ERROR('FINN height profile requires explicit MSL or AGL datum',RC)
       RETURN
    ENDIF
    ALLOCATE(Weights(HcoState%NX,HcoState%NY,HcoState%NZ), &
             Reason(HcoState%NX,HcoState%NY), SupportFraction(HcoState%NX,HcoState%NY), STAT=AS)
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
    REAL(hp) :: Support(HcoState%NX,HcoState%NY), Total(HcoState%NX,HcoState%NY)
    REAL(hp) :: MeanNum(HcoState%NX,HcoState%NY), TopNum(HcoState%NX,HcoState%NY)
    REAL(hp) :: M,T,G
    INTEGER :: I,J,S,Why
    Support=0.0_hp
    Total=0.0_hp
    MeanNum=0.0_hp
    TopNum=0.0_hp
    CALL HCO_EvalFld(HcoState,'FINNV25_GFAS_SUPPORT',Support,RC)
    IF (RC/=HCO_SUCCESS) RETURN
    CALL HCO_EvalFld(HcoState,'FINNV25_GFAS_TOTAL',Total,RC)
    IF (RC/=HCO_SUCCESS) RETURN
    CALL HCO_EvalFld(HcoState,'FINNV25_GFAS_MAMI_MOMENT',MeanNum,RC)
    IF (RC/=HCO_SUCCESS) RETURN
    CALL HCO_EvalFld(HcoState,'FINNV25_GFAS_APT_MOMENT',TopNum,RC)
    IF (RC/=HCO_SUCCESS) RETURN
    IF (ANY(.NOT.IEEE_IS_FINITE(Support)) .OR. ANY(.NOT.IEEE_IS_FINITE(Total)) .OR. &
        ANY(.NOT.IEEE_IS_FINITE(MeanNum)) .OR. ANY(.NOT.IEEE_IS_FINITE(TopNum))) THEN
       CALL HCO_ERROR('Nonfinite FINN GFAS height moment input',RC)
       RETURN
    ENDIF
    IF (ANY(Support<0.0_hp) .OR. ANY(Total<0.0_hp) .OR. &
        ANY(Support>Total*(1.0_hp+64.0_hp*EPSILON(1.0_hp)))) THEN
       CALL HCO_ERROR('Invalid FINN GFAS height support fraction',RC)
       RETURN
    ENDIF
    Reason=1
    Weights=0.0_hp
    SupportFraction=0.0_hp
    DO J=1,HcoState%NY
    DO I=1,HcoState%NX
       IF (Total(I,J)>0.0_hp) SupportFraction(I,J)=MIN(1.0_hp,Support(I,J)/Total(I,J))
       IF (Support(I,J)<=0.0_hp) CYCLE
       M=MeanNum(I,J)/Support(I,J)
       T=TopNum(I,J)/Support(I,J)
       G=0.0_hp
       IF (HeightDatum=='MSL') G=HcoState%Grid%ZSFC%Val(I,J)
       CALL FINN_HeightProfile(M,T,G,HcoState%Grid%BXHEIGHT_M%Val(I,J,:), &
                               .FALSE.,Weights(I,J,:),Why,S)
       IF (S==0) THEN
          Reason(I,J)=0
       ELSEIF (S==1 .AND. Why==2) THEN
          Reason(I,J)=2
       ELSEIF (S==3) THEN
          CALL HCO_ERROR('FINN GFAS requested height support exceeds model top',RC)
          RETURN
       ELSE
          CALL HCO_ERROR('Invalid FINN GFAS height ordering, datum, or model geometry',RC)
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
    A=Column*(1.0_hp-SupportFraction)
    CALL Diagn_Update(HcoState,cName='FINNProfileMissingSupport_'//TRIM(Name), &
                       AutoFill=0,Array2D=A,RC=RC)
    IF (RC/=HCO_SUCCESS) RETURN
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
    IF ( ALLOCATED(SupportFraction) ) DEALLOCATE(SupportFraction)
    HeightDatum='UNSET'
    IF ( ALLOCATED(Weights) ) DEALLOCATE(Weights)
    IF ( ALLOCATED(Reason) ) DEALLOCATE(Reason)
    ProfileEnabled=.FALSE.
    Fallback='pbl'
  END SUBROUTINE FINN_Profile_Final
END MODULE HCOX_FINN_PROFILE_MOD
