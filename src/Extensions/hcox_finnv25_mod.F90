!------------------------------------------------------------------------------
!                   Harmonized Emissions Component (HEMCO)
!------------------------------------------------------------------------------
! Apply the GFED vertical-injection method to species-resolved FINNv2.5
! fields named FINNV25_INJECT_<species> in HEMCO_Config.rc.
MODULE HCOX_FINNv25_MOD

  USE HCO_ERROR_MOD
  USE HCO_STATE_MOD,  ONLY : HCO_State
  USE HCOX_STATE_MOD, ONLY : Ext_State

  USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_FINITE
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: HCOX_FINNv25_Init
  PUBLIC :: HCOX_FINNv25_Run
  PUBLIC :: HCOX_FINNv25_Final
  PUBLIC :: HCOX_FINNv25_CloseOrigins

  INTEGER, SAVE                         :: ExtNrSaved = -1
  INTEGER, SAVE                         :: nSpc = 0
  INTEGER, ALLOCATABLE, SAVE            :: HcoIDs(:)
  CHARACTER(LEN=31), ALLOCATABLE, SAVE  :: SpcNames(:)
  REAL(hp), ALLOCATABLE, SAVE           :: SpcArr3D(:,:,:)
  REAL(hp), SAVE                        :: VerticalInjectFrac = 0.0_hp
  INTEGER, SAVE                         :: VerticalInjectLevels = 0
  LOGICAL, SAVE :: OriginTags = .FALSE.
  REAL(hp), ALLOCATABLE, SAVE :: OriginFire(:,:,:,:)
  CHARACTER(LEN=8), PARAMETER :: OriginParents(7) = [ CHARACTER(LEN=8) :: &
    'FSOAP','FSOAS','BRCSOA','NPBRCPOA','WTC','PBRCPOA','DBRCPOA' ]
  CHARACTER(LEN=3), PARAMETER :: OriginNames(3) = ['USA','CAN','ROW']

CONTAINS

  SUBROUTINE HCOX_FINNv25_Init( HcoState, ExtName, ExtState, RC )

    USE HCO_STATE_MOD,   ONLY : HCO_GetExtHcoID, HCO_GetHcoID
    USE HCO_EXTLIST_MOD, ONLY : GetExtNr, GetExtOpt
    USE HCOX_FIRE_INJECTION_MOD, ONLY : HCOX_FireInject_Validate

    TYPE(HCO_State), POINTER         :: HcoState
    CHARACTER(LEN=*), INTENT(IN)     :: ExtName
    TYPE(Ext_State), POINTER         :: ExtState
    INTEGER, INTENT(INOUT)           :: RC

    CHARACTER(LEN=255)               :: MSG, LOC
    INTEGER                          :: ExtNr, AS, S, O
    INTEGER, ALLOCATABLE             :: TmpIDs(:)
    CHARACTER(LEN=31), ALLOCATABLE   :: TmpNames(:)
    LOGICAL                          :: FOUND

    LOC = 'HCOX_FINNv25_Init (HCOX_FINNV25_MOD.F90)'
    ExtNr = GetExtNr( HcoState%Config%ExtList, TRIM(ExtName) )
    IF ( ExtNr <= 0 ) RETURN

    CALL HCO_ENTER( HcoState%Config%Err, LOC, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    IF ( ExtNrSaved > 0 ) THEN
       CALL HCO_ERROR( 'Only one FINNv2.5 injection instance is supported', &
                       RC, THISLOC=LOC )
       RETURN
    ENDIF

    CALL GetExtOpt( HcoState%Config, ExtNr,                            &
                    'FINNv25_vertical_injection_fraction',             &
                    OptValHp=VerticalInjectFrac, FOUND=FOUND, RC=RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    IF ( .NOT. FOUND ) VerticalInjectFrac = 0.0_hp

    CALL GetExtOpt( HcoState%Config, ExtNr,                            &
                    'FINNv25_vertical_injection_levels',               &
                    OptValInt=VerticalInjectLevels, FOUND=FOUND, RC=RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    IF ( .NOT. FOUND ) VerticalInjectLevels = 0

    CALL GetExtOpt(HcoState%Config,ExtNr,'FINNv25_BrC_origin_tags', &
                   OptValBool=OriginTags,FOUND=FOUND,RC=RC)
    IF (RC /= HCO_SUCCESS) RETURN
    IF (.NOT. FOUND) OriginTags=.FALSE.
    IF (OriginTags) THEN
      DO S=1,7
        DO O=1,3
          IF (HCO_GetHcoID(TRIM(OriginParents(S))//'_'//OriginNames(O),HcoState)<1) THEN
            CALL HCO_ERROR('FINNv25 origin tagging needs all 21 fire origin species',RC)
            RETURN
          ENDIF
        ENDDO
      ENDDO
      ALLOCATE(OriginFire(HcoState%NX,HcoState%NY,HcoState%NZ,7))
      OriginFire=0.0_hp
    ENDIF

    CALL HCOX_FireInject_Validate( VerticalInjectFrac,                    &
                                  VerticalInjectLevels, 'FINNv2.5', RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL HCO_GetExtHcoID( HcoState, ExtNr, TmpIDs, TmpNames, nSpc, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    IF ( nSpc < 1 ) THEN
       CALL HCO_ERROR( 'No FINNv2.5 injection species specified', RC, &
                       THISLOC=LOC )
       RETURN
    ENDIF

    ALLOCATE( HcoIDs(nSpc), SpcNames(nSpc),                         &
              SpcArr3D(HcoState%NX,HcoState%NY,HcoState%NZ), STAT=AS )
    IF ( AS /= 0 ) THEN
       CALL HCO_ERROR( 'Cannot allocate FINNv2.5 injection arrays', RC, &
                       THISLOC=LOC )
       RETURN
    ENDIF
    HcoIDs   = TmpIDs
    SpcNames = TmpNames
    DEALLOCATE( TmpIDs, TmpNames )

    ExtNrSaved       = ExtNr
    ExtState%FINNv25 = ExtNr
    IF ( VerticalInjectFrac > 0.0_hp ) THEN
       ExtState%FRAC_OF_PBL%DoUse = .TRUE.
    ENDIF

    IF ( HcoState%amIRoot ) THEN
       WRITE(MSG,*) 'FINNv2.5 injection species       : ', nSpc
       CALL HCO_MSG( MSG, LUN=HcoState%Config%hcoLogLUN )
       WRITE(MSG,*) 'FINNv2.5 vertical injection frac: ', VerticalInjectFrac
       CALL HCO_MSG( MSG, LUN=HcoState%Config%hcoLogLUN )
       WRITE(MSG,*) 'FINNv2.5 elevated levels        : ', VerticalInjectLevels
       CALL HCO_MSG( MSG, LUN=HcoState%Config%hcoLogLUN )
    ENDIF

    CALL HCO_LEAVE( HcoState%Config%Err, RC )
  END SUBROUTINE HCOX_FINNv25_Init


  SUBROUTINE HCOX_FINNv25_Run( ExtState, HcoState, RC )

    USE HCO_CALC_MOD,    ONLY : HCO_EvalFld
    USE HCO_FLUXARR_MOD, ONLY : HCO_EmisAdd
    USE HCOX_FIRE_INJECTION_MOD, ONLY : HCOX_FireInject_Apply
    USE HCO_STATE_MOD, ONLY : HCO_GetHcoID
    USE HCO_SCALE_MOD, ONLY : HCO_ScaleGet

    TYPE(Ext_State), POINTER         :: ExtState
    TYPE(HCO_State), POINTER         :: HcoState
    INTEGER, INTENT(INOUT)           :: RC

    CHARACTER(LEN=63)                :: FieldName
    CHARACTER(LEN=255)               :: MSG, LOC
    INTEGER :: N, S, O, TagId
    REAL(hp), TARGET                 :: SpcArr(HcoState%NX,HcoState%NY)
    REAL(hp), TARGET :: TagArr(HcoState%NX,HcoState%NY)
    REAL(hp) :: Parts(HcoState%NX,HcoState%NY,3)
    REAL(hp) :: InjectedSum(HcoState%NX,HcoState%NY,HcoState%NZ)
    REAL(hp) :: SourceErr, SourceScale

    LOC = 'HCOX_FINNv25_Run (HCOX_FINNV25_MOD.F90)'
    IF ( ExtState%FINNv25 <= 0 ) RETURN

    CALL HCO_ENTER( HcoState%Config%Err, LOC, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    IF (OriginTags) OriginFire=0.0_hp
    DO N = 1, nSpc
       IF ( HcoIDs(N) < 0 ) CYCLE
       FieldName = 'FINNV25_INJECT_' // TRIM(SpcNames(N))
       SpcArr = 0.0_hp
       CALL HCO_EvalFld( HcoState, TRIM(FieldName), SpcArr, RC )
       IF ( RC /= HCO_SUCCESS ) THEN
          MSG = 'Cannot evaluate FINNv2.5 field ' // TRIM(FieldName)
          CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
          RETURN
       ENDIF

       IF (OriginTags) THEN
          IF (ANY(.NOT. IEEE_IS_FINITE(SpcArr)) .OR. ANY(SpcArr<0.0_hp)) THEN
             CALL HCO_ERROR('Invalid FINNv25 parent origin source',RC)
             RETURN
          ENDIF
          IF (HcoState%Options%ScaleEmis .AND. HCO_ScaleGet(HcoIDs(N))/=1.0_hp) THEN
             CALL HCO_ERROR('Origin prototype requires unit universal emission scales',RC)
             RETURN
          ENDIF
       ENDIF
       S=0
       IF (OriginTags) THEN
          DO O=1,7
             IF (TRIM(SpcNames(N))==TRIM(OriginParents(O))) S=O
          ENDDO
          IF (S>0) THEN
             DO O=1,3
                FieldName='FINNV25_TAG_'//TRIM(SpcNames(N))//'_'//OriginNames(O)
                TagArr=0.0_hp
                CALL HCO_EvalFld(HcoState,TRIM(FieldName),TagArr,RC)
                IF (RC/=HCO_SUCCESS) RETURN
                IF (ANY(.NOT. IEEE_IS_FINITE(TagArr)) .OR. ANY(TagArr<0.0_hp)) THEN
                   CALL HCO_ERROR('Invalid native-partition FINNv25 origin source',RC)
                   RETURN
                ENDIF
                Parts(:,:,O)=TagArr
             ENDDO
             SourceErr=MAXVAL(ABS(SUM(Parts,DIM=3)-SpcArr))
             SourceScale=MAXVAL(ABS(SpcArr))
             IF (SourceErr>2.0e-6_hp*MAX(SourceScale,TINY(1.0_hp))) THEN
                CALL HCO_ERROR('FINNv25 country sources do not close to parent',RC)
                RETURN
             ENDIF
             IF (HcoState%amIRoot) WRITE(6,'(a,1x,a,2(1x,es22.14))') &
                'BRC_ORIGIN_SOURCE',TRIM(SpcNames(N)),SourceScale,SourceErr
          ENDIF
       ENDIF

       IF ( VerticalInjectFrac > 0.0_hp ) THEN
          CALL HCOX_FireInject_Apply( HcoState, ExtState, SpcArr,         &
                                     VerticalInjectFrac,                 &
                                     VerticalInjectLevels, SpcArr3D,     &
                                     'FINNv2.5', RC )
          IF ( RC /= HCO_SUCCESS ) RETURN

          CALL HCO_EmisAdd( HcoState, SpcArr3D, HcoIDs(N), RC,            &
                            ExtNr=ExtNrSaved )
       ELSE
          CALL HCO_EmisAdd( HcoState, SpcArr, HcoIDs(N), RC,              &
                            ExtNr=ExtNrSaved )
       ENDIF
       IF ( RC /= HCO_SUCCESS ) THEN
          MSG = 'HCO_EmisAdd error for FINNv2.5 ' // TRIM(SpcNames(N))
          CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
          RETURN
       ENDIF
       IF (OriginTags .AND. S>0) THEN
          InjectedSum=0.0_hp
          ! Cache the unchanged parent inventory flux. The untagged complement
          ! subtracts this actual inventory source, never a tag-sum correction.
          IF (VerticalInjectFrac>0.0_hp) THEN
             OriginFire(:,:,:,S)=SpcArr3D
          ELSE
             OriginFire(:,:,1,S)=SpcArr
          ENDIF
          DO O=1,3
             TagArr=Parts(:,:,O)
             TagId=HCO_GetHcoID(TRIM(SpcNames(N))//'_'//OriginNames(O),HcoState)
             IF (HcoState%Options%ScaleEmis .AND. HCO_ScaleGet(TagId)/=1.0_hp) THEN
                CALL HCO_ERROR('Origin prototype requires unit universal emission scales',RC)
                RETURN
             ENDIF
             IF (VerticalInjectFrac>0.0_hp) THEN
                CALL HCOX_FireInject_Apply(HcoState,ExtState,TagArr, &
                     VerticalInjectFrac,VerticalInjectLevels,SpcArr3D,'FINNv25 origin',RC)
                IF (RC/=HCO_SUCCESS) RETURN
                CALL HCO_EmisAdd(HcoState,SpcArr3D,TagId,RC,ExtNr=ExtNrSaved)
                InjectedSum=InjectedSum+SpcArr3D
             ELSE
                CALL HCO_EmisAdd(HcoState,TagArr,TagId,RC,ExtNr=ExtNrSaved)
                InjectedSum(:,:,1)=InjectedSum(:,:,1)+TagArr
             ENDIF
             IF (RC/=HCO_SUCCESS) RETURN
          ENDDO
          SourceErr=MAXVAL(ABS(InjectedSum-OriginFire(:,:,:,S)))
          SourceScale=MAXVAL(ABS(OriginFire(:,:,:,S)))
          IF (SourceErr>2.0e-6_hp*MAX(SourceScale,TINY(1.0_hp))) THEN
             CALL HCO_ERROR('FINNv25 injected origin sources do not close to parent',RC)
             RETURN
          ENDIF
          IF (HcoState%amIRoot) WRITE(6,'(a,1x,a,2(1x,es22.14))') &
             'BRC_ORIGIN_INJECTION',TRIM(SpcNames(N)),SourceScale,SourceErr
       ENDIF
    ENDDO

    CALL HCO_LEAVE( HcoState%Config%Err, RC )
  END SUBROUTINE HCOX_FINNv25_Run



  ! Called after ALL HEMCO extensions. Unrequested inventory/non-fire sources
  ! enter UNT. This never changes parents or geographically normalizes tags.
  SUBROUTINE HCOX_FINNv25_CloseOrigins(HcoState,RC)
    USE HCO_STATE_MOD, ONLY: HCO_GetHcoID
    USE HCO_FLUXARR_MOD, ONLY: HCO_EmisAdd
    USE HCO_SCALE_MOD, ONLY: HCO_ScaleGet
    TYPE(Hco_State), POINTER :: HcoState
    INTEGER, INTENT(OUT) :: RC
    INTEGER :: S,O,P,U,Count,Ids(7,4),ParentIds(7)
    REAL(hp) :: Flux(HcoState%NX,HcoState%NY,HcoState%NZ), Scale, Negative
    REAL(hp) :: OriginTotal(HcoState%NX,HcoState%NY,HcoState%NZ),Closure
    CHARACTER(LEN=3), PARAMETER :: AllOrigins(4)=['USA','CAN','ROW','UNT']
    RC=HCO_SUCCESS
    Count=0
    DO S=1,7
      ParentIds(S)=HCO_GetHcoID(TRIM(OriginParents(S)),HcoState)
      DO O=1,4
        Ids(S,O)=HCO_GetHcoID(TRIM(OriginParents(S))//'_'//AllOrigins(O),HcoState)
        IF (Ids(S,O)>0) Count=Count+1
      ENDDO
    ENDDO
    IF (Count==0) RETURN
    IF (Count/=28 .OR. ANY(ParentIds<1)) THEN
      CALL HCO_ERROR('Incomplete BrC origin species in HEMCO',RC)
      RETURN
    ENDIF
    DO S=1,7
      P=ParentIds(S); U=Ids(S,4); Flux=0.0_hp
      IF (HcoState%Options%ScaleEmis) THEN
        IF (HCO_ScaleGet(P)/=1.0_hp .OR. ANY([(HCO_ScaleGet(Ids(S,O))/=1.0_hp,O=1,4)])) THEN
          CALL HCO_ERROR('Origin prototype requires unit universal emission scales',RC)
          RETURN
        ENDIF
      ENDIF
      IF (ASSOCIATED(HcoState%Spc(P)%Emis)) THEN
        IF (ASSOCIATED(HcoState%Spc(P)%Emis%Val)) Flux=HcoState%Spc(P)%Emis%Val
      ENDIF
      IF (ANY(.NOT. IEEE_IS_FINITE(Flux))) THEN
        CALL HCO_ERROR('Nonfinite BrC parent emission flux',RC)
        RETURN
      ENDIF
      Scale=MAXVAL(ABS(Flux))
      IF (OriginTags) Flux=Flux-OriginFire(:,:,:,S)
      IF (ANY(.NOT. IEEE_IS_FINITE(Flux))) THEN
        CALL HCO_ERROR('Nonfinite BrC untagged emission flux',RC)
        RETURN
      ENDIF
      Negative=MAX(0.0_hp,-MINVAL(Flux))
      IF (Negative>2.0e-6_hp*MAX(Scale,TINY(1.0_hp))) THEN
        CALL HCO_ERROR('Requested FINNv25 fire source exceeds total parent emissions',RC)
        RETURN
      ENDIF
      ! Roundoff-only negative subtraction is reported, not used to reassign
      ! any geographic partition. UNT itself must have nonnegative emissions.
      IF (HcoState%amIRoot) WRITE(6,'(a,1x,a,2(1x,es22.14))') &
        'BRC_ORIGIN_UNTAGGED_SOURCE',TRIM(OriginParents(S)),MAXVAL(ABS(Flux)),Negative
      Flux=MAX(Flux,0.0_hp)
      CALL HCO_EmisAdd(HcoState,Flux,U,RC)
      IF (RC/=HCO_SUCCESS) RETURN
      ! Audit actual accumulated arrays after all source additions, including
      ! hp addition roundoff. Never adjust a geographic partition to close.
      OriginTotal=0.0_hp
      DO O=1,4
        IF (ASSOCIATED(HcoState%Spc(Ids(S,O))%Emis)) THEN
          IF (ASSOCIATED(HcoState%Spc(Ids(S,O))%Emis%Val)) &
            OriginTotal=OriginTotal+HcoState%Spc(Ids(S,O))%Emis%Val
        ENDIF
      ENDDO
      Flux=0.0_hp
      IF (ASSOCIATED(HcoState%Spc(P)%Emis)) THEN
        IF (ASSOCIATED(HcoState%Spc(P)%Emis%Val)) Flux=HcoState%Spc(P)%Emis%Val
      ENDIF
      IF (ANY(.NOT. IEEE_IS_FINITE(OriginTotal))) THEN
        CALL HCO_ERROR('Nonfinite accumulated BrC origin emissions',RC)
        RETURN
      ENDIF
      Closure=MAXVAL(ABS(OriginTotal-Flux))
      IF (HcoState%amIRoot) WRITE(6,'(a,1x,a,2(1x,es22.14))') &
        'BRC_ORIGIN_EMISSION_CLOSURE',TRIM(OriginParents(S)),Scale,Closure
      IF (Closure>2.0e-6_hp*MAX(Scale,TINY(1.0_hp))) THEN
        CALL HCO_ERROR('Accumulated BrC origin emissions do not close to parent',RC)
        RETURN
      ENDIF
    ENDDO
  END SUBROUTINE HCOX_FINNv25_CloseOrigins

  SUBROUTINE HCOX_FINNv25_Final( ExtState )
    TYPE(Ext_State), POINTER :: ExtState

    IF ( ALLOCATED(HcoIDs)   ) DEALLOCATE(HcoIDs)
    IF ( ALLOCATED(SpcNames) ) DEALLOCATE(SpcNames)
    IF (ALLOCATED(OriginFire)) DEALLOCATE(OriginFire)
    OriginTags=.FALSE.
    IF ( ALLOCATED(SpcArr3D) ) DEALLOCATE(SpcArr3D)
    ExtNrSaved           = -1
    nSpc                 = 0
    VerticalInjectFrac   = 0.0_hp
    VerticalInjectLevels = 0
    ExtState%FINNv25     = -1
  END SUBROUTINE HCOX_FINNv25_Final

END MODULE HCOX_FINNv25_MOD
