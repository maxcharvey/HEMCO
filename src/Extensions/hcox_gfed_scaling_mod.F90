!------------------------------------------------------------------------------
!                   Harmonized Emissions Component (HEMCO)                    !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hcox_gfed_scaling_mod.F90
!
! !DESCRIPTION: Module HCOX\_GFED\_SCALING\_MOD configures species derived
! from GFED CO so their final emitted mass-flux ratios are independent of the
! order of species in HEMCO\_Config.rc.
!
! !INTERFACE:
!
MODULE HCOX_GFED_SCALING_MOD
!
! !USES:
!
  USE HCO_PRECISION_MOD, ONLY : f4

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: CONFIGURE_GFED_CO_RATIOS
!
! !REVISION HISTORY:
!  06 Aug 2026 - M. Harvey - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
CONTAINS
!EOC
!------------------------------------------------------------------------------
!                   Harmonized Emissions Component (HEMCO)                    !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: CONFIGURE_GFED_CO_RATIOS
!
! !DESCRIPTION: Makes SOAP and FSOAP inherit the configured CO scalar and
! spatial scale field before their ratio is applied. Independent SOAP or
! FSOAP scaling is rejected because it would change the requested final ratio.
!
! !INTERFACE:
!
  SUBROUTINE CONFIGURE_GFED_CO_RATIOS( SpcNames, SOAPfrac, FSOAPfrac,       &
                                       SpcScal, SpcScaleField, NoScale,     &
                                       ErrMsg )
!
! !INPUT PARAMETERS:
!
    CHARACTER(LEN=*), INTENT(IN   ) :: SpcNames(:)
    REAL(f4),         INTENT(IN   ) :: SOAPfrac
    REAL(f4),         INTENT(IN   ) :: FSOAPfrac
    CHARACTER(LEN=*), INTENT(IN   ) :: NoScale
!
! !INPUT/OUTPUT PARAMETERS:
!
    REAL(f4),         INTENT(INOUT) :: SpcScal(:)
    CHARACTER(LEN=*), INTENT(INOUT) :: SpcScaleField(:)
!
! !OUTPUT PARAMETERS:
!
    CHARACTER(LEN=*), INTENT(  OUT) :: ErrMsg
!
! !REVISION HISTORY:
!  06 Aug 2026 - M. Harvey - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: CoIdx, FsoapIdx, N, SoapIdx

    ErrMsg  = ''
    CoIdx   = 0
    SoapIdx = 0
    FsoapIdx = 0

    IF ( SIZE(SpcNames) /= SIZE(SpcScal) .OR.                            &
         SIZE(SpcNames) /= SIZE(SpcScaleField) ) THEN
       ErrMsg = 'GFED CO-ratio scaling arrays have inconsistent sizes'
       RETURN
    ENDIF

    DO N = 1, SIZE(SpcNames)
       SELECT CASE ( TRIM(SpcNames(N)) )
          CASE ( 'CO' )
             IF ( CoIdx > 0 ) THEN
                ErrMsg = 'GFED CO-ratio scaling found duplicate CO species'
                RETURN
             ENDIF
             CoIdx = N
          CASE ( 'SOAP' )
             IF ( SoapIdx > 0 ) THEN
                ErrMsg = 'GFED CO-ratio scaling found duplicate SOAP species'
                RETURN
             ENDIF
             SoapIdx = N
          CASE ( 'FSOAP' )
             IF ( FsoapIdx > 0 ) THEN
                ErrMsg = 'GFED CO-ratio scaling found duplicate FSOAP species'
                RETURN
             ENDIF
             FsoapIdx = N
       END SELECT
    ENDDO

    CALL Validate_CO_Controls( 'SOAP', SoapIdx, SOAPfrac, CoIdx, SpcScal,  &
                               SpcScaleField, NoScale, ErrMsg )
    IF ( LEN_TRIM(ErrMsg) > 0 ) RETURN

    CALL Validate_CO_Controls( 'FSOAP', FsoapIdx, FSOAPfrac, CoIdx,        &
                               SpcScal, SpcScaleField, NoScale, ErrMsg )
    IF ( LEN_TRIM(ErrMsg) > 0 ) RETURN

    IF ( SoapIdx > 0 .AND. SOAPfrac > 0.0_f4 ) THEN
       SpcScal(SoapIdx)       = SpcScal(CoIdx)
       SpcScaleField(SoapIdx) = SpcScaleField(CoIdx)
    ENDIF
    IF ( FsoapIdx > 0 .AND. FSOAPfrac > 0.0_f4 ) THEN
       SpcScal(FsoapIdx)       = SpcScal(CoIdx)
       SpcScaleField(FsoapIdx) = SpcScaleField(CoIdx)
    ENDIF

  END SUBROUTINE CONFIGURE_GFED_CO_RATIOS
!EOC
!------------------------------------------------------------------------------
!                   Harmonized Emissions Component (HEMCO)                    !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Validate_CO_Controls
!
! !DESCRIPTION: Validates one CO-derived output before CO scaling is copied.
!
! !INTERFACE:
!
  SUBROUTINE Validate_CO_Controls( SpcName, SpcIdx, Fraction, CoIdx,       &
                                   SpcScal, SpcScaleField, NoScale, ErrMsg )
!
! !INPUT PARAMETERS:
!
    CHARACTER(LEN=*), INTENT(IN   ) :: SpcName
    INTEGER,          INTENT(IN   ) :: SpcIdx
    REAL(f4),         INTENT(IN   ) :: Fraction
    INTEGER,          INTENT(IN   ) :: CoIdx
    CHARACTER(LEN=*), INTENT(IN   ) :: NoScale
!
! !INPUT PARAMETERS:
!
    REAL(f4),         INTENT(IN   ) :: SpcScal(:)
    CHARACTER(LEN=*), INTENT(IN   ) :: SpcScaleField(:)
!
! !OUTPUT PARAMETERS:
!
    CHARACTER(LEN=*), INTENT(  OUT) :: ErrMsg
!
! !REVISION HISTORY:
!  06 Aug 2026 - M. Harvey - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC

    ErrMsg = ''
    IF ( SpcIdx < 1 .OR. Fraction <= 0.0_f4 ) RETURN

    IF ( CoIdx < 1 ) THEN
       ErrMsg = 'GFED CO-to-' // TRIM(SpcName) // ' ratio requires CO output'
       RETURN
    ENDIF

    IF ( SpcScal(SpcIdx) /= 1.0_f4 .AND.                              &
         SpcScal(SpcIdx) /= SpcScal(CoIdx) ) THEN
       ErrMsg = 'Scaling_' // TRIM(SpcName) //                           &
                ' conflicts with final GFED CO-to-' // TRIM(SpcName) //  &
                ' ratio'
       RETURN
    ENDIF

    IF ( TRIM(SpcScaleField(SpcIdx)) /= TRIM(NoScale) .AND.               &
         TRIM(SpcScaleField(SpcIdx)) /= TRIM(SpcScaleField(CoIdx)) ) THEN
       ErrMsg = 'ScaleField_' // TRIM(SpcName) //                        &
                ' conflicts with final GFED CO-to-' // TRIM(SpcName) //  &
                ' ratio'
       RETURN
    ENDIF

  END SUBROUTINE Validate_CO_Controls
!EOC
END MODULE HCOX_GFED_SCALING_MOD
