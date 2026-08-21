!------------------------------------------------------------------------------
!                   Harmonized Emissions Component (HEMCO)
!------------------------------------------------------------------------------
! Apply the GFED vertical-injection method to species-resolved FINNv2.5
! fields named FINNV25_INJECT_<species> in HEMCO_Config.rc.
MODULE HCOX_FINNv25_MOD

  USE HCO_ERROR_MOD
  USE HCO_STATE_MOD,  ONLY : HCO_State
  USE HCOX_STATE_MOD, ONLY : Ext_State

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: HCOX_FINNv25_Init
  PUBLIC :: HCOX_FINNv25_Run
  PUBLIC :: HCOX_FINNv25_Final

  INTEGER, SAVE                         :: ExtNrSaved = -1
  INTEGER, SAVE                         :: nSpc = 0
  INTEGER, ALLOCATABLE, SAVE            :: HcoIDs(:)
  CHARACTER(LEN=31), ALLOCATABLE, SAVE  :: SpcNames(:)
  REAL(hp), ALLOCATABLE, SAVE           :: SpcArr3D(:,:,:)
  REAL(hp), SAVE                        :: VerticalInjectFrac = 0.0_hp
  INTEGER, SAVE                         :: VerticalInjectLevels = 0

CONTAINS

  SUBROUTINE HCOX_FINNv25_Init( HcoState, ExtName, ExtState, RC )

    USE HCO_STATE_MOD,   ONLY : HCO_GetExtHcoID
    USE HCO_EXTLIST_MOD, ONLY : GetExtNr, GetExtOpt

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

    IF ( VerticalInjectFrac < 0.0_hp .OR. VerticalInjectFrac > 1.0_hp ) THEN
       CALL HCO_ERROR( 'FINNv2.5 injection fraction must be in [0,1]', &
                       RC, THISLOC=LOC )
       RETURN
    ENDIF
    IF ( VerticalInjectFrac == 0.0_hp .AND. VerticalInjectLevels /= 0 ) THEN
       CALL HCO_ERROR( 'FINNv2.5 zero injection requires zero levels', &
                       RC, THISLOC=LOC )
       RETURN
    ENDIF
    IF ( VerticalInjectFrac > 0.0_hp .AND. VerticalInjectLevels < 1 ) THEN
       CALL HCO_ERROR( 'FINNv2.5 positive injection requires levels', &
                       RC, THISLOC=LOC )
       RETURN
    ENDIF

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

    TYPE(Ext_State), POINTER         :: ExtState
    TYPE(HCO_State), POINTER         :: HcoState
    INTEGER, INTENT(INOUT)           :: RC

    CHARACTER(LEN=63)                :: FieldName
    CHARACTER(LEN=255)               :: MSG, LOC
    INTEGER                          :: I, J, L, N
    INTEGER                          :: PBL_MAX, N_FTLEV
    REAL(hp)                         :: F_OF_PBL, F_OF_FT, PBL_SUM
    REAL(hp)                         :: DELTPRES, TOTPRESFT
    REAL(hp), TARGET                 :: SpcArr(HcoState%NX,HcoState%NY)

    LOC = 'HCOX_FINNv25_Run (HCOX_FINNV25_MOD.F90)'
    IF ( ExtState%FINNv25 <= 0 ) RETURN

    CALL HCO_ENTER( HcoState%Config%Err, LOC, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

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

       IF ( VerticalInjectFrac > 0.0_hp ) THEN
          SpcArr3D = 0.0_hp
          DO J = 1, HcoState%NY
          DO I = 1, HcoState%NX
             IF ( SpcArr(I,J) <= 0.0_hp ) CYCLE

             PBL_MAX = 0
             DO L = HcoState%NZ, 1, -1
                IF ( ExtState%FRAC_OF_PBL%Arr%Val(I,J,L) > 0.0_hp ) THEN
                   PBL_MAX = L
                   EXIT
                ENDIF
             ENDDO

             PBL_SUM = 0.0_hp
             DO L = 1, PBL_MAX
                PBL_SUM = PBL_SUM + MAX( 0.0_hp,                         &
                                 ExtState%FRAC_OF_PBL%Arr%Val(I,J,L) )
             ENDDO
             IF ( PBL_SUM <= 0.0_hp ) THEN
                CALL HCO_ERROR( 'FINNv2.5 injection has invalid PBL fractions', &
                                RC, THISLOC=LOC )
                RETURN
             ENDIF

             N_FTLEV = VerticalInjectLevels
             IF ( PBL_MAX + N_FTLEV > HcoState%NZ ) THEN
                CALL HCO_ERROR( 'FINNv2.5 injection lacks requested FT levels', &
                                RC, THISLOC=LOC )
                RETURN
             ENDIF

             DO L = 1, PBL_MAX
                F_OF_PBL = MAX( 0.0_hp,                                  &
                         ExtState%FRAC_OF_PBL%Arr%Val(I,J,L) ) / PBL_SUM
                SpcArr3D(I,J,L) = SpcArr(I,J) *                           &
                           ( 1.0_hp - VerticalInjectFrac ) * F_OF_PBL
             ENDDO

             TOTPRESFT = HcoState%Grid%PEDGE%Val(I,J,PBL_MAX+1) -        &
                         HcoState%Grid%PEDGE%Val(I,J,PBL_MAX+N_FTLEV+1)
             IF ( TOTPRESFT <= 0.0_hp ) THEN
                CALL HCO_ERROR( 'FINNv2.5 injection has nonpositive FT depth', &
                                RC, THISLOC=LOC )
                RETURN
             ENDIF

             DO L = PBL_MAX+1, PBL_MAX+N_FTLEV
                DELTPRES = HcoState%Grid%PEDGE%Val(I,J,L) -               &
                           HcoState%Grid%PEDGE%Val(I,J,L+1)
                F_OF_FT = DELTPRES / TOTPRESFT
                SpcArr3D(I,J,L) = SpcArr(I,J) * VerticalInjectFrac * F_OF_FT
             ENDDO
          ENDDO
          ENDDO

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
