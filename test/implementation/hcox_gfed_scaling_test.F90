PROGRAM HCOX_GFED_SCALING_TEST

  USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY : IEEE_IS_FINITE, IEEE_QUIET_NAN, &
                                             IEEE_VALUE
  USE HCO_PRECISION_MOD,              ONLY : f4, hp
  USE HCOX_FIRE_INJECTION_MOD,        ONLY : HCOX_FireInject_Profile
  USE HCOX_GFED_SCALING_MOD,          ONLY : CONFIGURE_GFED_CO_RATIOS

  IMPLICIT NONE

  CHARACTER(LEN=31), PARAMETER :: BaseNames(4) =                        &
       (/ 'CO                             ',                            &
          'NAP                            ',                            &
          'SOAP                           ',                            &
          'FSOAP                          ' /)
  INTEGER :: Count, Order(4)
  LOGICAL :: Used(4)

  Count = 0
  Used  = .FALSE.
  CALL Permute( 1 )
  IF ( Count /= 24 ) ERROR STOP 1

  CALL Check_Isolation
  CALL Check_NAP_Isolation
  CALL Check_Mask_Propagation
  CALL Check_Repeated_Initialization
  CALL Check_Rejections
  CALL Check_Fire_Injection_Profile

  WRITE(*,'(a)') 'PASS: executable GFED CO-ratio scaling regression'

CONTAINS

  RECURSIVE SUBROUTINE Permute( Position )

    INTEGER, INTENT(IN) :: Position
    INTEGER :: Candidate

    IF ( Position > SIZE(Order) ) THEN
       CALL Check_Order
       Count = Count + 1
       RETURN
    ENDIF

    DO Candidate = 1, SIZE(Order)
       IF ( Used(Candidate) ) CYCLE
       Used(Candidate) = .TRUE.
       Order(Position) = Candidate
       CALL Permute( Position + 1 )
       Used(Candidate) = .FALSE.
    ENDDO

  END SUBROUTINE Permute

  SUBROUTINE Check_Order

    CHARACTER(LEN=255) :: ErrMsg
    CHARACTER(LEN=31)  :: Fields(4), Names(4)
    INTEGER            :: CoIdx, FsoapIdx, N, NapIdx, SoapIdx
    REAL(f4)           :: CoEmis, FsoapEmis, NapEmis, Scales(4), SoapEmis
    REAL(f4), PARAMETER :: Fraction = 0.013_f4
    ! Zero sentinel followed by the six production GFED4 CO fire-type factors.
    REAL(f4), PARAMETER :: RawCO(7) =                               &
         (/ 0.0_f4, 6.90e-2_f4, 1.21e-1_f4, 1.13e-1_f4,            &
            1.04e-1_f4, 2.60e-1_f4, 7.60e-2_f4 /)

    CoIdx    = 0
    FsoapIdx = 0
    NapIdx   = 0
    SoapIdx  = 0
    Fields   = 'none'
    DO N = 1, SIZE(Names)
       Names(N) = BaseNames(Order(N))
       SELECT CASE ( TRIM(Names(N)) )
          CASE ( 'CO' )
             Scales(N) = 1.05_f4
             CoIdx = N
          CASE ( 'NAP' )
             Scales(N) = 2.75e-4_f4
             NapIdx = N
          CASE ( 'SOAP' )
             Scales(N) = 1.0_f4
             SoapIdx = N
          CASE ( 'FSOAP' )
             Scales(N) = 1.0_f4
             FsoapIdx = N
       END SELECT
    ENDDO

    CALL CONFIGURE_GFED_CO_RATIOS( Names, Fraction, Fraction, Scales,     &
                                    Fields, 'none', ErrMsg )
    IF ( LEN_TRIM(ErrMsg) > 0 ) ERROR STOP 2

    IF ( Scales(NapIdx) /= 2.75e-4_f4 ) ERROR STOP 3
    IF ( Scales(SoapIdx) /= Scales(CoIdx) ) ERROR STOP 4
    IF ( Scales(FsoapIdx) /= Scales(CoIdx) ) ERROR STOP 5

    DO N = 1, SIZE(RawCO)
       CoEmis    = RawCO(N) * Scales(CoIdx)
       NapEmis   = RawCO(N) * Scales(NapIdx)
       SoapEmis  = RawCO(N) * Scales(SoapIdx) * Fraction
       FsoapEmis = RawCO(N) * Scales(FsoapIdx) * Fraction
       CALL Assert_Finite_Nonnegative( CoEmis, 6 )
       CALL Assert_Finite_Nonnegative( NapEmis, 7 )
       CALL Assert_Finite_Nonnegative( SoapEmis, 8 )
       CALL Assert_Finite_Nonnegative( FsoapEmis, 9 )
       CALL Assert_Close( CoEmis, RawCO(N) * 1.05_f4, 10 )
       CALL Assert_Close( NapEmis, RawCO(N) * 2.75e-4_f4, 11 )
       CALL Assert_Close( SoapEmis, RawCO(N) * 0.01365_f4, 12 )
       CALL Assert_Close( FsoapEmis, RawCO(N) * 0.01365_f4, 13 )
       IF ( CoEmis > 0.0_f4 ) THEN
          CALL Assert_Close( SoapEmis / CoEmis, Fraction, 14 )
          CALL Assert_Close( FsoapEmis / CoEmis, Fraction, 15 )
          CALL Assert_Close( NapEmis / CoEmis,                       &
                             2.75e-4_f4 / 1.05_f4, 16 )
       ELSE
          IF ( NapEmis /= 0.0_f4 .OR. SoapEmis /= 0.0_f4 .OR.       &
               FsoapEmis /= 0.0_f4 ) ERROR STOP 17
       ENDIF
    ENDDO

  END SUBROUTINE Check_Order

  SUBROUTINE Check_Isolation

    CHARACTER(LEN=255) :: ErrMsg
    CHARACTER(LEN=31)  :: Fields(4), Names(4)
    INTEGER            :: N
    REAL(f4)           :: CoScale, FsoapFraction, Scales(4), SoapFraction
    REAL(f4), PARAMETER :: CoScales(4) =                              &
         (/ 0.5_f4, 1.0_f4, 1.05_f4, 2.0_f4 /)

    Names         = BaseNames
    Fields        = 'none'
    SoapFraction  = 0.011_f4
    FsoapFraction = 0.017_f4

    DO N = 1, SIZE(CoScales)
       CoScale = CoScales(N)
       Scales  = (/ CoScale, 7.0e-4_f4, 1.0_f4, 1.0_f4 /)
       CALL CONFIGURE_GFED_CO_RATIOS( Names, SoapFraction, FsoapFraction, &
                                       Scales, Fields, 'none', ErrMsg )
       IF ( LEN_TRIM(ErrMsg) > 0 ) ERROR STOP 18
       CALL Assert_Close( Scales(1), CoScale, 19 )
       CALL Assert_Close( Scales(2), 7.0e-4_f4, 20 )
       CALL Assert_Close( Scales(3) * SoapFraction / Scales(1),      &
                          SoapFraction, 21 )
       CALL Assert_Close( Scales(4) * FsoapFraction / Scales(1),    &
                          FsoapFraction, 22 )
    ENDDO

    Names(1:3)  = (/ BaseNames(1), BaseNames(3), BaseNames(4) /)
    Fields(1:3) = 'none'
    Scales(1:3) = (/ 1.05_f4, 1.0_f4, 1.0_f4 /)
    CALL CONFIGURE_GFED_CO_RATIOS( Names(1:3), SoapFraction,            &
                                    FsoapFraction, Scales(1:3),          &
                                    Fields(1:3), 'none', ErrMsg )
    IF ( LEN_TRIM(ErrMsg) > 0 ) ERROR STOP 23
    CALL Assert_Close( Scales(2), Scales(1), 24 )
    CALL Assert_Close( Scales(3), Scales(1), 25 )

  END SUBROUTINE Check_Isolation

  SUBROUTINE Check_NAP_Isolation

    CHARACTER(LEN=255) :: ErrMsg
    CHARACTER(LEN=31)  :: Fields(4), Names(4)
    INTEGER            :: N
    REAL(f4)           :: CoEmis, FsoapEmis, NapEmis, Scales(4), SoapEmis
    REAL(f4), PARAMETER :: NapScales(4) =                             &
         (/ 0.0_f4, 2.75e-4_f4, 7.0e-4_f4, 2.0_f4 /)

    Names  = BaseNames
    Fields = 'none'
    DO N = 1, SIZE(NapScales)
       Scales = (/ 1.05_f4, NapScales(N), 1.0_f4, 1.0_f4 /)
       CALL CONFIGURE_GFED_CO_RATIOS( Names, 0.013_f4, 0.013_f4, Scales, &
                                       Fields, 'none', ErrMsg )
       IF ( LEN_TRIM(ErrMsg) > 0 ) ERROR STOP 26
       CoEmis    = 3.0_f4 * Scales(1)
       NapEmis   = 3.0_f4 * Scales(2)
       SoapEmis  = 3.0_f4 * Scales(3) * 0.013_f4
       FsoapEmis = 3.0_f4 * Scales(4) * 0.013_f4
       CALL Assert_Close( CoEmis, 3.15_f4, 27 )
       CALL Assert_Close( NapEmis, 3.0_f4 * NapScales(N), 28 )
       CALL Assert_Close( SoapEmis, 0.04095_f4, 29 )
       CALL Assert_Close( FsoapEmis, 0.04095_f4, 30 )
    ENDDO

  END SUBROUTINE Check_NAP_Isolation

  SUBROUTINE Check_Mask_Propagation

    CHARACTER(LEN=255) :: ErrMsg
    CHARACTER(LEN=31)  :: Fields(4), Names(4)
    REAL(f4)           :: Scales(4)

    Names     = BaseNames
    Fields    = 'none'
    Fields(1) = 'CO_MASK'
    Scales    = (/ 1.05_f4, 2.75e-4_f4, 1.0_f4, 1.0_f4 /)
    CALL CONFIGURE_GFED_CO_RATIOS( Names, 0.013_f4, 0.013_f4, Scales,     &
                                    Fields, 'none', ErrMsg )
    IF ( LEN_TRIM(ErrMsg) > 0 ) ERROR STOP 31
    IF ( TRIM(Fields(2)) /= 'none' ) ERROR STOP 32
    IF ( TRIM(Fields(3)) /= 'CO_MASK' ) ERROR STOP 33
    IF ( TRIM(Fields(4)) /= 'CO_MASK' ) ERROR STOP 34
    CALL CONFIGURE_GFED_CO_RATIOS( Names, 0.013_f4, 0.013_f4, Scales,     &
                                    Fields, 'none', ErrMsg )
    IF ( LEN_TRIM(ErrMsg) > 0 ) ERROR STOP 35

  END SUBROUTINE Check_Mask_Propagation

  SUBROUTINE Check_Repeated_Initialization

    CHARACTER(LEN=255) :: ErrMsg
    CHARACTER(LEN=31)  :: Fields(4), Names(4)
    REAL(f4)           :: Scales(4)

    Names  = BaseNames
    Fields = 'none'
    Scales = (/ 1.05_f4, 2.75e-4_f4, 1.0_f4, 1.0_f4 /)
    CALL CONFIGURE_GFED_CO_RATIOS( Names, 0.013_f4, 0.013_f4, Scales,     &
                                    Fields, 'none', ErrMsg )
    IF ( LEN_TRIM(ErrMsg) > 0 ) ERROR STOP 36
    CALL CONFIGURE_GFED_CO_RATIOS( Names, 0.013_f4, 0.013_f4, Scales,     &
                                    Fields, 'none', ErrMsg )
    IF ( LEN_TRIM(ErrMsg) > 0 ) ERROR STOP 37
    CALL Assert_Close( Scales(3), 1.05_f4, 38 )
    CALL Assert_Close( Scales(4), 1.05_f4, 39 )

  END SUBROUTINE Check_Repeated_Initialization

  SUBROUTINE Check_Fire_Injection_Profile

    INTEGER :: Status
    REAL(hp) :: Flux, Frac, PBL(20), PEdge(21), Profile(20)

    PBL = 0.0_hp
    PBL(1:3) = (/ 0.2_hp, 0.3_hp, 0.5_hp /)
    DO Status = 1, SIZE(PEdge)
       PEdge(Status) = 1000.0_hp - 25.0_hp * REAL(Status-1,hp)
    ENDDO

    Flux = 10.0_hp
    Frac = 0.35_hp
    CALL HCOX_FireInject_Profile( Flux, PBL, PEdge, Frac, 10, Profile, Status )
    IF ( Status /= 0 ) ERROR STOP 60
    CALL Assert_Close_Hp( SUM(Profile), Flux, 61 )
    CALL Assert_Close_Hp( SUM(Profile(1:3)), Flux*(1.0_hp-Frac), 62 )
    CALL Assert_Close_Hp( SUM(Profile(4:13)), Flux*Frac, 63 )

    CALL HCOX_FireInject_Profile( Flux, PBL, PEdge, Frac, 15, Profile, Status )
    IF ( Status /= 0 ) ERROR STOP 64
    CALL Assert_Close_Hp( SUM(Profile), Flux, 65 )
    CALL Assert_Close_Hp( SUM(Profile(4:18)), Flux*Frac, 66 )

    PBL = 0.0_hp
    PBL(19) = 1.0_hp
    CALL HCOX_FireInject_Profile( Flux, PBL, PEdge, Frac, 2, Profile, Status )
    IF ( Status /= 2 ) ERROR STOP 67
    CALL Assert_Close_Hp( SUM(Profile), 0.0_hp, 68 )

    PBL = 0.0_hp
    CALL HCOX_FireInject_Profile( 0.0_hp, PBL, PEdge, Frac, 10, Profile, Status )
    IF ( Status /= 0 ) ERROR STOP 69
    CALL Assert_Close_Hp( SUM(Profile), 0.0_hp, 70 )

    CALL HCOX_FireInject_Profile( Flux, PBL, PEdge, Frac, 10, Profile, Status )
    IF ( Status /= 1 ) ERROR STOP 71
    CALL Assert_Close_Hp( SUM(Profile), 0.0_hp, 72 )

    PBL(1) = 1.0_hp
    CALL HCOX_FireInject_Profile( Flux, PBL, PEdge, 0.0_hp, 0, Profile, Status )
    IF ( Status /= 0 ) ERROR STOP 73
    CALL Assert_Close_Hp( Profile(1), Flux, 74 )
    CALL Assert_Close_Hp( SUM(Profile(2:)), 0.0_hp, 75 )
    CALL HCOX_FireInject_Profile( Flux, PBL, PEdge,                         &
                                   IEEE_VALUE(Flux, IEEE_QUIET_NAN), 10,     &
                                   Profile, Status )
    IF ( Status /= 4 ) ERROR STOP 76
    CALL HCOX_FireInject_Profile( -Flux, PBL, PEdge, Frac, 10, Profile, Status )
    IF ( Status /= 4 ) ERROR STOP 77
    CALL HCOX_FireInject_Profile( Flux, PBL, PEdge, 0.0_hp, 10, Profile, Status )
    IF ( Status /= 4 ) ERROR STOP 78
    CALL HCOX_FireInject_Profile( Flux, PBL, PEdge, Frac, -1, Profile, Status )
    IF ( Status /= 4 ) ERROR STOP 79

  END SUBROUTINE Check_Fire_Injection_Profile

  SUBROUTINE Assert_Close_Hp( Actual, Expected, Code )

    REAL(hp), INTENT(IN) :: Actual, Expected
    INTEGER, INTENT(IN)  :: Code

    IF ( ABS(Actual-Expected) > 1.0e-12_hp ) ERROR STOP Code

  END SUBROUTINE Assert_Close_Hp

  SUBROUTINE Check_Rejections

    CHARACTER(LEN=255) :: ErrMsg
    CHARACTER(LEN=31)  :: Fields(4), FieldsBefore(4), Names(4)
    REAL(f4)           :: Scales(4), ScalesBefore(4)

    Names  = BaseNames
    Fields = 'none'
    Scales = (/ 1.05_f4, 2.75e-4_f4, 2.0_f4, 1.0_f4 /)
    FieldsBefore = Fields
    ScalesBefore = Scales
    CALL CONFIGURE_GFED_CO_RATIOS( Names, 0.013_f4, 0.013_f4, Scales,     &
                                    Fields, 'none', ErrMsg )
    IF ( INDEX(ErrMsg, 'Scaling_SOAP') < 1 ) ERROR STOP 40
    CALL Assert_Unchanged( Scales, ScalesBefore, Fields, FieldsBefore, 41 )

    Fields = 'none'
    Scales = (/ 1.05_f4, 2.75e-4_f4, 1.0_f4, 2.0_f4 /)
    FieldsBefore = Fields
    ScalesBefore = Scales
    CALL CONFIGURE_GFED_CO_RATIOS( Names, 0.013_f4, 0.013_f4, Scales,     &
                                    Fields, 'none', ErrMsg )
    IF ( INDEX(ErrMsg, 'Scaling_FSOAP') < 1 ) ERROR STOP 42
    CALL Assert_Unchanged( Scales, ScalesBefore, Fields, FieldsBefore, 43 )

    Fields    = 'none'
    Fields(4) = 'FSOAP_MASK'
    Scales    = (/ 1.05_f4, 2.75e-4_f4, 1.0_f4, 1.0_f4 /)
    FieldsBefore = Fields
    ScalesBefore = Scales
    CALL CONFIGURE_GFED_CO_RATIOS( Names, 0.013_f4, 0.013_f4, Scales,     &
                                    Fields, 'none', ErrMsg )
    IF ( INDEX(ErrMsg, 'ScaleField_FSOAP') < 1 ) ERROR STOP 44
    CALL Assert_Unchanged( Scales, ScalesBefore, Fields, FieldsBefore, 45 )

    Names(1:3)  = (/ BaseNames(2), BaseNames(3), BaseNames(4) /)
    Fields(1:3) = 'none'
    Scales(1:3) = (/ 2.75e-4_f4, 1.0_f4, 1.0_f4 /)
    FieldsBefore(1:3) = Fields(1:3)
    ScalesBefore(1:3) = Scales(1:3)
    CALL CONFIGURE_GFED_CO_RATIOS( Names(1:3), 0.013_f4, 0.013_f4,        &
                                    Scales(1:3), Fields(1:3), 'none', ErrMsg )
    IF ( INDEX(ErrMsg, 'requires CO output') < 1 ) ERROR STOP 46
    CALL Assert_Unchanged( Scales(1:3), ScalesBefore(1:3), Fields(1:3), &
                           FieldsBefore(1:3), 47 )

    Names  = (/ BaseNames(1), BaseNames(1), BaseNames(3), BaseNames(4) /)
    Fields = 'none'
    Scales = (/ 1.05_f4, 1.05_f4, 1.0_f4, 1.0_f4 /)
    FieldsBefore = Fields
    ScalesBefore = Scales
    CALL CONFIGURE_GFED_CO_RATIOS( Names, 0.013_f4, 0.013_f4, Scales,     &
                                    Fields, 'none', ErrMsg )
    IF ( INDEX(ErrMsg, 'duplicate CO') < 1 ) ERROR STOP 48
    CALL Assert_Unchanged( Scales, ScalesBefore, Fields, FieldsBefore, 49 )

    Names  = BaseNames
    Fields = 'none'
    Scales = (/ 1.05_f4, 2.75e-4_f4, 1.0_f4, 1.0_f4 /)
    FieldsBefore = Fields
    ScalesBefore = Scales
    CALL CONFIGURE_GFED_CO_RATIOS( Names, 0.013_f4, 0.013_f4,             &
                                    Scales(1:3), Fields, 'none', ErrMsg )
    IF ( INDEX(ErrMsg, 'inconsistent sizes') < 1 ) ERROR STOP 50
    CALL Assert_Unchanged( Scales, ScalesBefore, Fields, FieldsBefore, 51 )

  END SUBROUTINE Check_Rejections

  SUBROUTINE Assert_Unchanged( Scales, ExpectedScales, Fields,             &
                               ExpectedFields, Code )

    REAL(f4),         INTENT(IN) :: Scales(:), ExpectedScales(:)
    CHARACTER(LEN=*), INTENT(IN) :: Fields(:), ExpectedFields(:)
    INTEGER,          INTENT(IN) :: Code

    IF ( ANY(Scales /= ExpectedScales) .OR. ANY(Fields /= ExpectedFields) ) THEN
       WRITE(*,'(a,i0)') 'FAIL mutation on rejected configuration ', Code
       ERROR STOP 97
    ENDIF

  END SUBROUTINE Assert_Unchanged

  SUBROUTINE Assert_Close( Actual, Expected, Code )

    REAL(f4), INTENT(IN) :: Actual, Expected
    INTEGER,  INTENT(IN) :: Code
    REAL(f4) :: Tolerance

    Tolerance = MAX( 1.0e-20_f4, 1.0e-6_f4 * ABS(Expected) )
    IF ( ABS(Actual - Expected) > Tolerance ) THEN
       WRITE(*,'(a,i0,2(a,es15.7))') 'FAIL assertion ', Code,             &
            ': actual=', Actual, ' expected=', Expected
       ERROR STOP 98
    ENDIF

  END SUBROUTINE Assert_Close

  SUBROUTINE Assert_Finite_Nonnegative( Value, Code )

    REAL(f4), INTENT(IN) :: Value
    INTEGER,  INTENT(IN) :: Code

    IF ( .NOT. IEEE_IS_FINITE(Value) .OR. Value < 0.0_f4 ) THEN
       WRITE(*,'(a,i0,a,es15.7)') 'FAIL finite/nonnegative assertion ',   &
            Code, ': value=', Value
       ERROR STOP 99
    ENDIF

  END SUBROUTINE Assert_Finite_Nonnegative

END PROGRAM HCOX_GFED_SCALING_TEST
