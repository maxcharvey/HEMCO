!------------------------------------------------------------------------------
!                   Harmonized Emissions Component (HEMCO)
!------------------------------------------------------------------------------
! Read species-resolved GFAS fields and apply the common fire-injection profile.
MODULE HCOX_GFAS_MOD

  USE HCO_ERROR_MOD
  USE HCO_STATE_MOD,  ONLY : HCO_State
  USE HCOX_STATE_MOD, ONLY : Ext_State

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: HCOX_GFAS_Init
  PUBLIC :: HCOX_GFAS_Run
  PUBLIC :: HCOX_GFAS_Final

  INTEGER, SAVE                         :: ExtNrSaved = -1
  INTEGER, SAVE                         :: nSpc = 0
  INTEGER, ALLOCATABLE, SAVE            :: HcoIDs(:)
  CHARACTER(LEN=31), ALLOCATABLE, SAVE  :: SpcNames(:)
  REAL(hp), ALLOCATABLE, SAVE           :: SpcArr3D(:,:,:)
  REAL(hp), SAVE                        :: VerticalInjectFrac = 0.0_hp
  INTEGER, SAVE                         :: VerticalInjectLevels = 0

CONTAINS

  SUBROUTINE HCOX_GFAS_Init( HcoState, ExtName, ExtState, RC )

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

    LOC = 'HCOX_GFAS_Init (HCOX_GFAS_MOD.F90)'
    ExtNr = GetExtNr( HcoState%Config%ExtList, TRIM(ExtName) )
    IF ( ExtNr <= 0 ) RETURN

    CALL HCO_ENTER( HcoState%Config%Err, LOC, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    IF ( ExtNrSaved > 0 ) THEN
       CALL HCO_ERROR( 'Only one GFAS injection instance is supported', &
                       RC, THISLOC=LOC )
       RETURN
    ENDIF

    CALL GetExtOpt( HcoState%Config, ExtNr,                            &
                    'GFAS_vertical_injection_fraction',                &
                    OptValHp=VerticalInjectFrac, FOUND=FOUND, RC=RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    IF ( .NOT. FOUND ) VerticalInjectFrac = 0.0_hp

    CALL GetExtOpt( HcoState%Config, ExtNr,                            &
                    'GFAS_vertical_injection_levels',                  &
                    OptValInt=VerticalInjectLevels, FOUND=FOUND, RC=RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    IF ( .NOT. FOUND ) VerticalInjectLevels = 0

    CALL HCOX_FireInject_Validate( VerticalInjectFrac,                 &
                                  VerticalInjectLevels, 'GFAS', RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL HCO_GetExtHcoID( HcoState, ExtNr, TmpIDs, TmpNames, nSpc, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    IF ( nSpc < 1 ) THEN
       CALL HCO_ERROR( 'No GFAS injection species specified', RC, THISLOC=LOC )
       RETURN
    ENDIF

    ALLOCATE( HcoIDs(nSpc), SpcNames(nSpc),                         &
              SpcArr3D(HcoState%NX,HcoState%NY,HcoState%NZ), STAT=AS )
    IF ( AS /= 0 ) THEN
       CALL HCO_ERROR( 'Cannot allocate GFAS injection arrays', RC, THISLOC=LOC )
       RETURN
    ENDIF
    HcoIDs   = TmpIDs
    SpcNames = TmpNames
    DEALLOCATE( TmpIDs, TmpNames )

    ExtNrSaved       = ExtNr
    ExtState%GFASInject = ExtNr
    IF ( VerticalInjectFrac > 0.0_hp ) THEN
       ExtState%FRAC_OF_PBL%DoUse = .TRUE.
    ENDIF

    IF ( HcoState%amIRoot ) THEN
       WRITE(MSG,*) 'GFAS injection species       : ', nSpc
       CALL HCO_MSG( MSG, LUN=HcoState%Config%hcoLogLUN )
       WRITE(MSG,*) 'GFAS vertical injection frac: ', VerticalInjectFrac
       CALL HCO_MSG( MSG, LUN=HcoState%Config%hcoLogLUN )
       WRITE(MSG,*) 'GFAS elevated levels        : ', VerticalInjectLevels
       CALL HCO_MSG( MSG, LUN=HcoState%Config%hcoLogLUN )
    ENDIF

    CALL HCO_LEAVE( HcoState%Config%Err, RC )
  END SUBROUTINE HCOX_GFAS_Init


  SUBROUTINE HCOX_GFAS_Run( ExtState, HcoState, RC )

    USE HCO_CALC_MOD,    ONLY : HCO_EvalFld
    USE HCO_FLUXARR_MOD, ONLY : HCO_EmisAdd
    USE HCOX_FIRE_INJECTION_MOD, ONLY : HCOX_FireInject_Apply

    TYPE(Ext_State), POINTER         :: ExtState
    TYPE(HCO_State), POINTER         :: HcoState
    INTEGER, INTENT(INOUT)           :: RC

    CHARACTER(LEN=63)                :: FieldName
    CHARACTER(LEN=255)               :: MSG, LOC
    INTEGER                          :: N
    REAL(hp), TARGET                 :: SpcArr(HcoState%NX,HcoState%NY)
    REAL(hp), TARGET                 :: SpcArr2(HcoState%NX,HcoState%NY)

    LOC = 'HCOX_GFAS_Run (HCOX_GFAS_MOD.F90)'
    IF ( ExtState%GFASInject <= 0 ) RETURN

    CALL HCO_ENTER( HcoState%Config%Err, LOC, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    DO N = 1, nSpc
       IF ( HcoIDs(N) < 0 ) CYCLE
       SpcArr = 0.0_hp
       IF ( TRIM(SpcNames(N)) == 'PRPE' ) THEN
          CALL HCO_EvalFld( HcoState, 'GFAS_INJECT_PRPE1', SpcArr, RC )
          IF ( RC /= HCO_SUCCESS ) RETURN
          SpcArr2 = 0.0_hp
          CALL HCO_EvalFld( HcoState, 'GFAS_INJECT_PRPE2', SpcArr2, RC )
          IF ( RC /= HCO_SUCCESS ) RETURN
          SpcArr = SpcArr + SpcArr2
       ELSE
          FieldName = 'GFAS_INJECT_' // TRIM(SpcNames(N))
          CALL HCO_EvalFld( HcoState, TRIM(FieldName), SpcArr, RC )
          IF ( RC /= HCO_SUCCESS ) THEN
             MSG = 'Cannot evaluate GFAS field ' // TRIM(FieldName)
             CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
             RETURN
          ENDIF
       ENDIF

       CALL HCOX_FireInject_Apply( HcoState, ExtState, SpcArr,            &
                                  VerticalInjectFrac,                    &
                                  VerticalInjectLevels, SpcArr3D,        &
                                  'GFAS', RC )
       IF ( RC /= HCO_SUCCESS ) RETURN
       CALL HCO_EmisAdd( HcoState, SpcArr3D, HcoIDs(N), RC,              &
                         ExtNr=ExtNrSaved )
       IF ( RC /= HCO_SUCCESS ) THEN
          MSG = 'HCO_EmisAdd error for GFAS ' // TRIM(SpcNames(N))
          CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
          RETURN
       ENDIF
    ENDDO

    CALL HCO_LEAVE( HcoState%Config%Err, RC )
  END SUBROUTINE HCOX_GFAS_Run


  SUBROUTINE HCOX_GFAS_Final( ExtState )
    TYPE(Ext_State), POINTER :: ExtState

    IF ( ALLOCATED(HcoIDs)   ) DEALLOCATE(HcoIDs)
    IF ( ALLOCATED(SpcNames) ) DEALLOCATE(SpcNames)
    IF ( ALLOCATED(SpcArr3D) ) DEALLOCATE(SpcArr3D)
    ExtNrSaved           = -1
    nSpc                 = 0
    VerticalInjectFrac   = 0.0_hp
    VerticalInjectLevels = 0
    ExtState%GFASInject  = -1
  END SUBROUTINE HCOX_GFAS_Final

END MODULE HCOX_GFAS_MOD
