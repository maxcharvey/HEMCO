! Pure, conservative kernels for opt-in FINNv2.5 auxiliary fire profiles.
MODULE HCOX_FINN_PROFILE_KERNEL_MOD
  USE HCO_PRECISION_MOD, ONLY: hp => f8
  USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_FINITE
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: FINN_Normalize, FINN_Allocate
CONTAINS
  PURE SUBROUTINE FINN_Normalize( Input, Weight, Status )
    REAL(hp), INTENT(IN) :: Input(:)
    REAL(hp), INTENT(OUT) :: Weight(:)
    INTEGER, INTENT(OUT) :: Status
    REAL(hp) :: Peak, Total
    INTEGER :: K
    Weight = 0.0_hp
    Status = 2
    IF ( SIZE(Input) /= SIZE(Weight) .OR. SIZE(Input) < 1 ) RETURN
    IF ( ANY(.NOT. IEEE_IS_FINITE(Input)) .OR. ANY(Input < 0.0_hp) ) RETURN
    Status = 1 ! recognized empty support, caller chooses explicit fallback
    Peak = MAXVAL(Input)
    IF ( Peak <= 0.0_hp ) RETURN
    Weight = Input / Peak
    Total = SUM(Weight)
    Weight = Weight / Total
    K = MAXLOC(Weight, DIM=1)
    Weight(K) = 0.0_hp
    Weight(K) = 1.0_hp - SUM(Weight)
    Status = 0
  END SUBROUTINE FINN_Normalize

  PURE SUBROUTINE FINN_Allocate( Flux, Weight, Profile, Status )
    REAL(hp), INTENT(IN) :: Flux, Weight(:)
    REAL(hp), INTENT(OUT) :: Profile(:)
    INTEGER, INTENT(OUT) :: Status
    INTEGER :: K
    Profile = 0.0_hp
    Status = 2
    IF ( SIZE(Profile) /= SIZE(Weight) .OR. SIZE(Weight) < 1 ) RETURN
    IF ( .NOT. IEEE_IS_FINITE(Flux) .OR. Flux < 0.0_hp ) RETURN
    IF ( ANY(.NOT. IEEE_IS_FINITE(Weight)) .OR. ANY(Weight < 0.0_hp) ) RETURN
    IF ( ABS(SUM(Weight)-1.0_hp) > 64.0_hp*EPSILON(1.0_hp) ) RETURN
    Profile = Flux * Weight
    ! Place only rounding residue in the largest bin, never a physical tail.
    K = MAXLOC(Weight, DIM=1)
    Profile(K) = 0.0_hp
    Profile(K) = Flux - SUM(Profile)
    IF ( ANY(Profile < 0.0_hp) .OR. ANY(.NOT. IEEE_IS_FINITE(Profile)) ) RETURN
    Status = 0
  END SUBROUTINE FINN_Allocate
END MODULE HCOX_FINN_PROFILE_KERNEL_MOD
