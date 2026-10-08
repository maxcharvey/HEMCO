! Source allocation only: never repairs transported origin inventories.
MODULE HCOX_ORIGIN_SOURCE_KERNEL_MOD
  USE HCO_PRECISION_MOD, ONLY: fp
  USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_FINITE
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: ORIGIN_SOURCE_ALLOCATE
CONTAINS
  PURE SUBROUTINE ORIGIN_SOURCE_ALLOCATE(Parent,Raw,Injected,Parts,Status,RelativeError)
    REAL(fp), INTENT(IN) :: Parent,Raw(3),Injected
    REAL(fp), INTENT(INOUT) :: Parts(3)
    INTEGER, INTENT(OUT) :: Status
    REAL(fp), INTENT(OUT) :: RelativeError
    REAL(fp) :: Total,Trial(3),Tolerance,TrialSum
    LOGICAL :: SafeSum
    Status=1;RelativeError=0.0_fp
    IF (.NOT. IEEE_IS_FINITE(Parent) .OR. .NOT. IEEE_IS_FINITE(Injected)) RETURN
    IF (ANY(.NOT. IEEE_IS_FINITE(Raw))) RETURN
    IF (Parent<0.0_fp .OR. Injected<0.0_fp .OR. ANY(Raw<0.0_fp)) RETURN
    CALL NONNEG_SUM3(Raw,Total,SafeSum)
    IF (.NOT. SafeSum) RETURN
    IF (Parent==0.0_fp) THEN
      IF (Total>0.0_fp .OR. Injected>0.0_fp) THEN
        Status=2
        RETURN
      ENDIF
      Parts=0.0_fp;Status=0
      RETURN
    ENDIF
    IF (Total==0.0_fp) THEN
      Status=2
      RETURN
    ENDIF
    RelativeError=ABS(Total-Parent)/MAX(Total,Parent)
    IF (RelativeError>2.0e-6_fp) THEN
      Status=3
      RETURN
    ENDIF
    ! Direct products preserve rare labels that subtraction would erase.
    Trial=Injected*(Raw/Total)
    IF (ANY(.NOT. IEEE_IS_FINITE(Trial)) .OR. ANY(Trial<0.0_fp)) RETURN
    ! Unsupported labels and representational loss both refuse explicitly.
    IF (ANY((Raw==0.0_fp) .AND. (Trial/=0.0_fp))) THEN
      Status=4
      RETURN
    ENDIF
    IF (Injected>0.0_fp .AND. ANY((Raw>0.0_fp) .AND. (Trial==0.0_fp))) THEN
      Status=4
      RETURN
    ENDIF
    Tolerance=5.0e-5_fp
    IF (STORAGE_SIZE(Parent)>32) Tolerance=5.0e-12_fp
    IF (Injected>0.0_fp) THEN
      CALL NONNEG_SUM3(Trial,TrialSum,SafeSum)
      IF (.NOT. SafeSum) THEN
        Status=5
        RETURN
      ENDIF
      IF (ABS(TrialSum-Injected)/Injected>Tolerance) THEN
        Status=5
        RETURN
      ENDIF
    ENDIF
    Parts=Trial;Status=0
  END SUBROUTINE ORIGIN_SOURCE_ALLOCATE
  PURE SUBROUTINE NONNEG_SUM3(Values,Total,Valid)
    REAL(fp), INTENT(IN) :: Values(3)
    REAL(fp), INTENT(OUT) :: Total
    LOGICAL, INTENT(OUT) :: Valid
    INTEGER :: O
    Valid=.FALSE.;Total=0.0_fp
    DO O=1,3
      IF (Values(O)>HUGE(1.0_fp)-Total) RETURN
      Total=Total+Values(O)
    ENDDO
    Valid=IEEE_IS_FINITE(Total)
  END SUBROUTINE NONNEG_SUM3
END MODULE HCOX_ORIGIN_SOURCE_KERNEL_MOD
