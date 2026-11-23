! Pure, conservative kernels for opt-in FINNv2.5 auxiliary fire profiles.
MODULE HCOX_FINN_PROFILE_KERNEL_MOD
  USE HCO_PRECISION_MOD, ONLY: hp => f8
  USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_FINITE
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: FINN_Normalize, FINN_Allocate, FINN_Scale_Height
CONTAINS
  ! Conservatively scale an already normalized model-grid profile in metres
  ! above the current model surface. Each input bin is uniform per unit height;
  ! its mass moves from [zlo,zhi] to Scale*[zlo,zhi]. Layer HEIGHTS are depths,
  ! not MSL altitudes. Never subtract terrain or multiply emission mass.
  ! Status: 0 success, 2 invalid input, 3 positive mass above model top.
  ! Overflow is the fraction outside the domain. Overflow fails closed: there
  ! is no silent clipping, accumulation at the top, or tail renormalization.
  PURE SUBROUTINE FINN_Scale_Height( Input, Heights, Scale, Weight, Overflow, Status )
    REAL(hp), INTENT(IN) :: Input(:), Heights(:), Scale
    REAL(hp), INTENT(OUT) :: Weight(:), Overflow
    INTEGER, INTENT(OUT) :: Status
    REAL(hp) :: Edge(0:SIZE(Input)), Bottom, Top, Width, Overlap
    REAL(hp) :: Total, Correction
    INTEGER :: N, K, L, First, Peak
    Weight = 0.0_hp
    Overflow = 0.0_hp
    Status = 2
    N = SIZE(Input)
    IF ( N < 1 .OR. SIZE(Weight) /= N .OR. SIZE(Heights) /= N ) RETURN
    IF ( .NOT. IEEE_IS_FINITE(Scale) ) RETURN
    IF ( Scale <= 0.0_hp ) RETURN
    IF ( ANY(.NOT. IEEE_IS_FINITE(Input)) .OR. ANY(Input < 0.0_hp) ) RETURN
    IF ( ANY(Input > 1.0_hp+64.0_hp*EPSILON(1.0_hp)) ) RETURN
    Total = SUM(Input)
    IF ( ABS(Total-1.0_hp) > 64.0_hp*EPSILON(1.0_hp) ) RETURN
    ! Identity must not smooth the profile or change rounding. Geometry is
    ! intentionally unused for this case, as in the production bypass.
    IF ( Scale == 1.0_hp ) THEN
       Weight = Input
       Status = 0
       RETURN
    ENDIF
    IF ( ANY(.NOT. IEEE_IS_FINITE(Heights)) .OR. ANY(Heights <= 0.0_hp) ) RETURN
    Edge(0) = 0.0_hp
    DO K = 1, N
       IF ( Heights(K) > HUGE(1.0_hp)-Edge(K-1) ) RETURN
       Edge(K) = Edge(K-1)+Heights(K)
       IF ( Edge(K) <= Edge(K-1) ) RETURN
    ENDDO
    IF ( Scale > 1.0_hp ) THEN
       IF ( Edge(N) > HUGE(1.0_hp)/Scale ) RETURN
    ENDIF
    First = 1
    DO K = 1, N
       IF ( Input(K) == 0.0_hp ) CYCLE
       Bottom = Scale*Edge(K-1)
       Top = Scale*Edge(K)
       Width = Top-Bottom
       IF ( Width <= 0.0_hp ) THEN
          Weight = 0.0_hp
          RETURN
       ENDIF
       IF ( Top > Edge(N) ) THEN
          Overflow = Overflow+Input(K)*MIN(1.0_hp,(Top-MAX(Bottom,Edge(N)))/Width)
       ENDIF
       DO WHILE ( First < N )
          IF ( Edge(First) > Bottom ) EXIT
          First = First+1
       ENDDO
       DO L = First, N
          IF ( Edge(L-1) >= Top ) EXIT
          Overlap = MAX(0.0_hp,MIN(Top,Edge(L))-MAX(Bottom,Edge(L-1)))
          Weight(L) = Weight(L)+Input(K)*(Overlap/Width)
       ENDDO
    ENDDO
    IF ( Overflow > 0.0_hp ) THEN
       Weight = 0.0_hp
       Status = 3
       RETURN
    ENDIF
    ! Only assign floating-point residue; do not repair missing physical mass.
    Correction = Total-SUM(Weight)
    IF ( ABS(Correction) > 128.0_hp*EPSILON(1.0_hp) ) THEN
       Weight = 0.0_hp
       RETURN
    ENDIF
    Peak = MAXLOC(Weight,DIM=1)
    Weight(Peak) = Weight(Peak)+Correction
    IF ( ANY(Weight < 0.0_hp) .OR. ANY(.NOT. IEEE_IS_FINITE(Weight)) ) THEN
       Weight = 0.0_hp
       RETURN
    ENDIF
    Status = 0
  END SUBROUTINE FINN_Scale_Height

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
