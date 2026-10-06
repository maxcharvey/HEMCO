!------------------------------------------------------------------------------
!                   Harmonized Emissions Component (HEMCO)
!------------------------------------------------------------------------------
! Apply the GFED vertical-injection method to species-resolved FINNv2.5
! fields named FINNV25_INJECT_<species> in HEMCO_Config.rc.
MODULE HCOX_FINNv25_MOD

  USE HCO_ERROR_MOD
  USE HCO_STATE_MOD,  ONLY : HCO_State
  USE HCOX_STATE_MOD, ONLY : Ext_State

  USE HCOX_FINN_PROFILE_MOD

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: HCOX_FINNv25_Init
  PUBLIC :: HCOX_FINNv25_Run
  PUBLIC :: HCOX_FINNv25_Final

  INTEGER, SAVE                         :: ExtNrSaved = -1
  INTEGER, SAVE                         :: nSpc = 0
  LOGICAL, SAVE                         :: InputsChecked = .FALSE.
  INTEGER, ALLOCATABLE, SAVE            :: HcoIDs(:)
  CHARACTER(LEN=31), ALLOCATABLE, SAVE  :: SpcNames(:)
  REAL(hp), ALLOCATABLE, SAVE           :: SpcArr3D(:,:,:)
  REAL(hp), SAVE                        :: VerticalInjectFrac = 0.0_hp
  INTEGER, SAVE                         :: VerticalInjectLevels = 0

CONTAINS

  SUBROUTINE HCOX_FINNv25_Init( HcoState, ExtName, ExtState, RC )

    USE HCO_STATE_MOD,   ONLY : HCO_GetExtHcoID
    USE HCO_EXTLIST_MOD, ONLY : GetExtNr, GetExtOpt
    USE HCOX_FIRE_INJECTION_MOD, ONLY : HCOX_FireInject_Validate

    TYPE(HCO_State), POINTER         :: HcoState
    CHARACTER(LEN=*), INTENT(IN)     :: ExtName
    TYPE(Ext_State), POINTER         :: ExtState
    INTEGER, INTENT(INOUT)           :: RC

    CHARACTER(LEN=255)               :: MSG, LOC
    INTEGER                          :: ExtNr, AS
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

    CALL FINN_Profile_Init(HcoState,ExtState,ExtNr,HcoIDs,SpcNames,RC)
    IF ( RC /= HCO_SUCCESS ) RETURN
    IF ( ProfileEnabled .AND. (VerticalInjectFrac /= 0.0_hp .OR. VerticalInjectLevels /= 0) ) THEN
       CALL HCO_ERROR('FINN GFAS profile requires legacy fraction/levels both zero', RC)
       RETURN
    ENDIF

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
    USE HCO_DATACONT_MOD, ONLY : ListCont_Find
    USE HCO_TYPES_MOD, ONLY : ListCont, HCO_UFLAG_ONCE
    USE HCO_FLUXARR_MOD, ONLY : HCO_EmisAdd
    USE HCOX_FIRE_INJECTION_MOD, ONLY : HCOX_FireInject_Apply

    TYPE(Ext_State), POINTER         :: ExtState
    TYPE(HCO_State), POINTER         :: HcoState
    INTEGER, INTENT(INOUT)           :: RC

    CHARACTER(LEN=63)                :: FieldName
    CHARACTER(LEN=255)               :: MSG, LOC
    INTEGER                          :: N
    TYPE(ListCont), POINTER           :: InputCont
    LOGICAL                          :: Found
    REAL(hp), TARGET                 :: SpcArr(HcoState%NX,HcoState%NY)

    LOC = 'HCOX_FINNv25_Run (HCOX_FINNV25_MOD.F90)'
    IF ( ExtState%FINNv25 <= 0 ) RETURN

    CALL HCO_ENTER( HcoState%Config%Err, LOC, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    IF ( .NOT. InputsChecked ) THEN
       ! Validate the effective containers, not only our generated templates.
       ! A daily annual file with EF/EY remains on the initialization slice.
       DO N=1,nSpc+1
          IF ( N <= nSpc ) THEN
             IF ( HcoIDs(N) < 0 ) CYCLE
             FieldName = 'FINNV25_INJECT_'//TRIM(SpcNames(N))
          ELSE
             IF ( .NOT. ProfileEnabled ) CYCLE
             FieldName = 'FINNV25_GFAS_REFERENCE'
          ENDIF
          InputCont => NULL()
          CALL ListCont_Find(HcoState%EmisList,TRIM(FieldName),Found,InputCont)
          IF ( .NOT. Found ) THEN
             CALL HCO_ERROR('Missing FINNv2.5 input container '//TRIM(FieldName),RC)
             RETURN
          ENDIF
          IF ( InputCont%Dct%Dta%ncRead .AND. &
               InputCont%Dct%Dta%UpdtFlag == HCO_UFLAG_ONCE ) THEN
             CALL HCO_ERROR('FINNv2.5 daily input uses read-once time flag: '// &
                            TRIM(FieldName)//'; use RF for FINN or EFY for GFAS',RC)
             RETURN
          ENDIF
       ENDDO
       InputsChecked = .TRUE.
    ENDIF

    IF ( ProfileEnabled ) THEN
       CALL FINN_Profile_Prepare(HcoState,RC)
       IF ( RC /= HCO_SUCCESS ) RETURN
    ENDIF

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

       IF ( ProfileEnabled ) THEN
          CALL FINN_Profile_Apply(HcoState,ExtState,SpcArr,SpcArr3D,RC)
          IF ( RC /= HCO_SUCCESS ) RETURN
          CALL HCO_EmisAdd(HcoState,SpcArr3D,HcoIDs(N),RC,ExtNr=ExtNrSaved)
          IF ( RC /= HCO_SUCCESS ) RETURN
          CALL FINN_Profile_Report(HcoState,ExtState,HcoIDs(N),SpcNames(N),SpcArr,SpcArr3D,RC)
       ELSEIF ( VerticalInjectFrac > 0.0_hp ) THEN
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
    ENDDO

    CALL HCO_LEAVE( HcoState%Config%Err, RC )
  END SUBROUTINE HCOX_FINNv25_Run


  SUBROUTINE HCOX_FINNv25_Final( ExtState )
    TYPE(Ext_State), POINTER :: ExtState

    CALL FINN_Profile_Final()
    InputsChecked = .FALSE.
    IF ( ALLOCATED(HcoIDs)   ) DEALLOCATE(HcoIDs)
    IF ( ALLOCATED(SpcNames) ) DEALLOCATE(SpcNames)
    IF ( ALLOCATED(SpcArr3D) ) DEALLOCATE(SpcArr3D)
    ExtNrSaved           = -1
    nSpc                 = 0
    VerticalInjectFrac   = 0.0_hp
    VerticalInjectLevels = 0
    ExtState%FINNv25     = -1
  END SUBROUTINE HCOX_FINNv25_Final

END MODULE HCOX_FINNv25_MOD
