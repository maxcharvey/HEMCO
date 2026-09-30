! Pure, conservative kernels for opt-in FINNv2.5 auxiliary fire profiles.
MODULE HCOX_FINN_PROFILE_KERNEL_MOD
  USE HCO_PRECISION_MOD, ONLY: hp => f8
  USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_FINITE
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: FINN_Normalize, FINN_Allocate, FINN_HeightProfile
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

  PURE SUBROUTINE FINN_HeightProfile( MeanMSL, TopMSL, Ground, Depth, Gaussian, Weight, Why, Status )
    REAL(hp), INTENT(IN) :: MeanMSL,TopMSL,Ground,Depth(:)
    LOGICAL, INTENT(IN) :: Gaussian
    REAL(hp), INTENT(OUT) :: Weight(:)
    INTEGER, INTENT(OUT) :: Why,Status
    REAL(hp) :: Mu,Top,H,Sigma,Lo,Hi,A,B,Edge(SIZE(Depth)+1),Raw(SIZE(Depth))
    INTEGER :: K
    Weight=0.0_hp
    Why=0
    Status=2
    IF (SIZE(Depth)<1 .OR. SIZE(Weight)/=SIZE(Depth)) RETURN
    IF (.NOT.IEEE_IS_FINITE(MeanMSL) .OR. .NOT.IEEE_IS_FINITE(TopMSL) .OR. &
        .NOT.IEEE_IS_FINITE(Ground) .OR. ANY(.NOT.IEEE_IS_FINITE(Depth))) RETURN
    IF (ANY(Depth<=0.0_hp) .OR. TopMSL<MeanMSL) RETURN
    Edge(1)=0.0_hp
    DO K=1,SIZE(Depth)
       Edge(K+1)=Edge(K)+Depth(K)
    ENDDO
    IF (ANY(.NOT.IEEE_IS_FINITE(Edge))) RETURN
    Mu=MeanMSL-Ground
    Top=TopMSL-Ground
    IF (.NOT.IEEE_IS_FINITE(Mu) .OR. .NOT.IEEE_IS_FINITE(Top)) RETURN
    IF (Mu<=0.0_hp) THEN
       Why=2 ! terrain-incompatible height; caller applies explicit fallback
       Status=1
       RETURN
    ENDIF
    Lo=0.0_hp
    Hi=Mu
    H=0.0_hp
    IF (Gaussian) THEN
       H=MIN(Top-Mu,Mu)
       Lo=Mu-H
       Hi=Mu+H
    ENDIF
    IF (Hi>Edge(SIZE(Edge))) THEN
       Status=3 ! no silent top truncation or top-layer remainder dumping
       RETURN
    ENDIF
    IF (Gaussian .AND. H==0.0_hp) THEN
       DO K=1,SIZE(Depth)
          ! Point exactly on an interior edge belongs to the upper layer.
          IF (Mu<Edge(K+1) .OR. K==SIZE(Depth)) THEN
             Weight(K)=1.0_hp
             Status=0
             RETURN
          ENDIF
       ENDDO
    ENDIF
    Raw=0.0_hp
    Sigma=H/3.0_hp
    DO K=1,SIZE(Depth)
       A=MAX(Edge(K),Lo)
       B=MIN(Edge(K+1),Hi)
       IF (B<=A) CYCLE
       IF (Gaussian) THEN
          Raw(K)=0.5_hp*(ERF((B-Mu)/(SQRT(2.0_hp)*Sigma))- &
                         ERF((A-Mu)/(SQRT(2.0_hp)*Sigma)))
       ELSE
          Raw(K)=B-A
       ENDIF
    ENDDO
    CALL FINN_Normalize(Raw,Weight,Status)
  END SUBROUTINE FINN_HeightProfile
END MODULE HCOX_FINN_PROFILE_KERNEL_MOD
