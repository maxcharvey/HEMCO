!------------------------------------------------------------------------------
!                   Harmonized Emissions Component (HEMCO)
!------------------------------------------------------------------------------
! Shared, mass-conserving vertical allocation for biomass-burning inventories.
MODULE HCOX_FIRE_INJECTION_MOD

  USE HCO_ERROR_MOD
  USE HCO_STATE_MOD,  ONLY : HCO_State
  USE HCOX_STATE_MOD, ONLY : Ext_State

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: HCOX_FireInject_Validate
  PUBLIC :: HCOX_FireInject_Apply

CONTAINS

  SUBROUTINE HCOX_FireInject_Validate( ElevatedFrac, ElevatedLevels, &
                                       InventoryName, RC )

    REAL(hp), INTENT(IN)         :: ElevatedFrac
    INTEGER, INTENT(IN)          :: ElevatedLevels
    CHARACTER(LEN=*), INTENT(IN) :: InventoryName
    INTEGER, INTENT(INOUT)       :: RC

    CHARACTER(LEN=255)           :: MSG

    IF ( ElevatedFrac < 0.0_hp .OR. ElevatedFrac > 1.0_hp ) THEN
       MSG = TRIM(InventoryName) // ' elevated fraction must be in [0,1]'
       CALL HCO_ERROR( MSG, RC )
       RETURN
    ENDIF
    IF ( ElevatedFrac == 0.0_hp .AND. ElevatedLevels /= 0 ) THEN
       MSG = TRIM(InventoryName) // ' zero elevated fraction requires zero levels'
       CALL HCO_ERROR( MSG, RC )
       RETURN
    ENDIF
    IF ( ElevatedFrac > 0.0_hp .AND. ElevatedLevels < 1 ) THEN
       MSG = TRIM(InventoryName) // ' positive elevated fraction requires levels'
       CALL HCO_ERROR( MSG, RC )
       RETURN
    ENDIF

  END SUBROUTINE HCOX_FireInject_Validate


  SUBROUTINE HCOX_FireInject_Apply( HcoState, ExtState, Flux2D,      &
                                    ElevatedFrac, ElevatedLevels,    &
                                    Flux3D, InventoryName, RC )

    TYPE(HCO_State), POINTER       :: HcoState
    TYPE(Ext_State), POINTER       :: ExtState
    REAL(hp), INTENT(IN)           :: Flux2D(:,:)
    REAL(hp), INTENT(IN)           :: ElevatedFrac
    INTEGER, INTENT(IN)            :: ElevatedLevels
    REAL(hp), INTENT(OUT)          :: Flux3D(:,:,:)
    CHARACTER(LEN=*), INTENT(IN)   :: InventoryName
    INTEGER, INTENT(INOUT)         :: RC

    CHARACTER(LEN=255)             :: MSG
    INTEGER                        :: I, J, L, PBL_MAX
    REAL(hp)                       :: DELTPRES, F_OF_FT, F_OF_PBL
    REAL(hp)                       :: PBL_SUM, TOTPRESFT

    Flux3D = 0.0_hp
    DO J = 1, HcoState%NY
    DO I = 1, HcoState%NX
       IF ( Flux2D(I,J) <= 0.0_hp ) CYCLE

       PBL_MAX = 0
       DO L = HcoState%NZ, 1, -1
          IF ( ExtState%FRAC_OF_PBL%Arr%Val(I,J,L) > 0.0_hp ) THEN
             PBL_MAX = L
             EXIT
          ENDIF
       ENDDO

       PBL_SUM = 0.0_hp
       DO L = 1, PBL_MAX
          PBL_SUM = PBL_SUM + MAX( 0.0_hp,                              &
                                   ExtState%FRAC_OF_PBL%Arr%Val(I,J,L) )
       ENDDO
       IF ( PBL_SUM <= 0.0_hp ) THEN
          MSG = TRIM(InventoryName) // ' injection has invalid PBL fractions'
          CALL HCO_ERROR( MSG, RC )
          RETURN
       ENDIF

       IF ( PBL_MAX + ElevatedLevels > HcoState%NZ ) THEN
          MSG = TRIM(InventoryName) // ' injection lacks requested FT levels'
          CALL HCO_ERROR( MSG, RC )
          RETURN
       ENDIF

       DO L = 1, PBL_MAX
          F_OF_PBL = MAX( 0.0_hp,                                       &
                          ExtState%FRAC_OF_PBL%Arr%Val(I,J,L) ) / PBL_SUM
          Flux3D(I,J,L) = Flux2D(I,J) * ( 1.0_hp - ElevatedFrac ) * F_OF_PBL
       ENDDO

       TOTPRESFT = HcoState%Grid%PEDGE%Val(I,J,PBL_MAX+1) -             &
                   HcoState%Grid%PEDGE%Val(I,J,PBL_MAX+ElevatedLevels+1)
       IF ( TOTPRESFT <= 0.0_hp ) THEN
          MSG = TRIM(InventoryName) // ' injection has nonpositive FT depth'
          CALL HCO_ERROR( MSG, RC )
          RETURN
       ENDIF

       DO L = PBL_MAX+1, PBL_MAX+ElevatedLevels
          DELTPRES = HcoState%Grid%PEDGE%Val(I,J,L) -                    &
                     HcoState%Grid%PEDGE%Val(I,J,L+1)
          F_OF_FT = DELTPRES / TOTPRESFT
          Flux3D(I,J,L) = Flux2D(I,J) * ElevatedFrac * F_OF_FT
       ENDDO
    ENDDO
    ENDDO

  END SUBROUTINE HCOX_FireInject_Apply

END MODULE HCOX_FIRE_INJECTION_MOD
