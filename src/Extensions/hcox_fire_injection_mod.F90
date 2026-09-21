!------------------------------------------------------------------------------
!                   Harmonized Emissions Component (HEMCO)
!------------------------------------------------------------------------------
! Shared, mass-conserving vertical allocation for biomass-burning inventories.
MODULE HCOX_FIRE_INJECTION_MOD

  USE HCO_ERROR_MOD
  USE HCO_STATE_MOD,  ONLY : HCO_State
  USE HCOX_STATE_MOD, ONLY : Ext_State
  USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY : IEEE_IS_FINITE

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: HCOX_FireInject_Validate
  PUBLIC :: HCOX_FireInject_Apply
  PUBLIC :: HCOX_FireInject_Profile

CONTAINS

  SUBROUTINE HCOX_FireInject_Validate( ElevatedFrac, ElevatedLevels, &
                                       InventoryName, RC )

    REAL(hp), INTENT(IN)         :: ElevatedFrac
    INTEGER, INTENT(IN)          :: ElevatedLevels
    CHARACTER(LEN=*), INTENT(IN) :: InventoryName
    INTEGER, INTENT(INOUT)       :: RC

    CHARACTER(LEN=255)           :: MSG

    ! The inclusive comparison also rejects NaN, for which both ordinary
    ! ordered comparisons are false.
    IF ( .NOT. ( ElevatedFrac >= 0.0_hp .AND. &
                 ElevatedFrac <= 1.0_hp ) ) THEN
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
    INTEGER                        :: I, J, Status
    REAL(hp)                       :: PBLFrac(HcoState%NZ)
    REAL(hp)                       :: PEdge(HcoState%NZ+1)
    REAL(hp)                       :: Profile(HcoState%NZ)

    Flux3D = 0.0_hp

    ! The legacy zero-control path is surface-only and must not require the
    ! PBL-fraction field. Validate source fluxes before any field access.
    IF ( ANY( .NOT. IEEE_IS_FINITE(Flux2D) ) .OR. &
         ANY( Flux2D < 0.0_hp ) ) THEN
       MSG = TRIM(InventoryName) // ' injection has invalid source flux'
       CALL HCO_ERROR( MSG, RC )
       RETURN
    ENDIF
    IF ( ElevatedFrac == 0.0_hp .AND. ElevatedLevels == 0 ) THEN
       Flux3D(:,:,1) = Flux2D
       RETURN
    ENDIF

    DO J = 1, HcoState%NY
    DO I = 1, HcoState%NX
       IF ( Flux2D(I,J) == 0.0_hp ) CYCLE

       PBLFrac = ExtState%FRAC_OF_PBL%Arr%Val(I,J,:)
       PEdge   = HcoState%Grid%PEDGE%Val(I,J,:)
       CALL HCOX_FireInject_Profile( Flux2D(I,J), PBLFrac, PEdge,        &
                                     ElevatedFrac, ElevatedLevels,       &
                                     Profile, Status )
       IF ( Status /= 0 ) THEN
          SELECT CASE ( Status )
             CASE ( 1 )
                MSG = TRIM(InventoryName) // ' injection has invalid PBL fractions'
             CASE ( 2 )
                MSG = TRIM(InventoryName) // ' injection lacks requested FT levels'
             CASE ( 3 )
                MSG = TRIM(InventoryName) // ' injection has nonpositive FT depth'
             CASE DEFAULT
                MSG = TRIM(InventoryName) // ' injection has invalid numeric input'
          END SELECT
          CALL HCO_ERROR( MSG, RC )
          RETURN
       ENDIF
       Flux3D(I,J,:) = Profile
    ENDDO
    ENDDO

  END SUBROUTINE HCOX_FireInject_Apply


  PURE SUBROUTINE HCOX_FireInject_Profile( Flux, PBLFrac, PEdge,          &
                                           ElevatedFrac, ElevatedLevels,  &
                                           Profile, Status )

    REAL(hp), INTENT(IN)  :: Flux, PBLFrac(:), PEdge(:), ElevatedFrac
    INTEGER, INTENT(IN)   :: ElevatedLevels
    REAL(hp), INTENT(OUT) :: Profile(:)
    INTEGER, INTENT(OUT)  :: Status

    INTEGER  :: L, PBLMax
    REAL(hp) :: DeltaPres, PBLTotal, FTTotal

    Profile = 0.0_hp
    Status  = 0
    IF ( SIZE(Profile) /= SIZE(PBLFrac) .OR. &
         SIZE(PEdge) /= SIZE(PBLFrac) + 1 ) THEN
       Status = 1
       RETURN
    ENDIF
    IF ( .NOT. IEEE_IS_FINITE(Flux) .OR. &
         .NOT. IEEE_IS_FINITE(ElevatedFrac) .OR. &
         ANY( .NOT. IEEE_IS_FINITE(PBLFrac) ) .OR. &
         ANY( .NOT. IEEE_IS_FINITE(PEdge) ) ) THEN
       Status = 4
       RETURN
    ENDIF
    IF ( Flux == 0.0_hp ) RETURN
    IF ( ElevatedFrac == 0.0_hp .AND. ElevatedLevels == 0 ) THEN
       Profile(1) = Flux
       RETURN
    ENDIF
    IF ( Flux < 0.0_hp .OR. ElevatedLevels < 1 .OR. &
         ElevatedFrac < 0.0_hp .OR. &
         ElevatedFrac > 1.0_hp .OR. ANY( PBLFrac < 0.0_hp ) .OR. &
         ANY( PEdge(:SIZE(PEdge)-1) <= PEdge(2:) ) ) THEN
       Status = 4
       RETURN
   ENDIF

    PBLMax = 0
    DO L = SIZE(PBLFrac), 1, -1
       IF ( PBLFrac(L) > 0.0_hp ) THEN
          PBLMax = L
          EXIT
       ENDIF
    ENDDO
    PBLTotal = SUM( MAX( 0.0_hp, PBLFrac(1:PBLMax) ) )
    IF ( PBLTotal <= 0.0_hp ) THEN
       Status = 1
       RETURN
    ENDIF
    IF ( PBLMax + ElevatedLevels > SIZE(Profile) ) THEN
       Status = 2
       RETURN
    ENDIF

    Profile(1:PBLMax) = Flux * ( 1.0_hp - ElevatedFrac ) *               &
                          MAX( 0.0_hp, PBLFrac(1:PBLMax) ) / PBLTotal
    FTTotal = PEdge(PBLMax+1) - PEdge(PBLMax+ElevatedLevels+1)
    IF ( FTTotal <= 0.0_hp ) THEN
       Status = 3
       Profile = 0.0_hp
       RETURN
    ENDIF
    DO L = PBLMax+1, PBLMax+ElevatedLevels
       DeltaPres = PEdge(L) - PEdge(L+1)
       Profile(L) = Flux * ElevatedFrac * DeltaPres / FTTotal
    ENDDO

  END SUBROUTINE HCOX_FireInject_Profile

END MODULE HCOX_FIRE_INJECTION_MOD
