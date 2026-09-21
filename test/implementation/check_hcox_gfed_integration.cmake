file(READ "${GFED_FILE}" GFED_SOURCE)
file(READ "${GFED_INCLUDE_FILE}" GFED_INCLUDE_SOURCE)
file(READ "${SCALING_FILE}" SCALING_SOURCE)

set(REQUIRED_GFED_TEXT
    "CALL CONFIGURE_GFED_CO_RATIOS("
    "IF ( TRIM(SpcName) == 'SOAP' ) SpcName = 'CO'"
    "IF ( TRIM(SpcName) == 'FSOAP'   ) SpcName = 'CO'"
    "CALL HCOX_FireInject_Apply( HcoState, ExtState, SpcArr,"
    "CALL HCO_EmisAdd( HcoState, Inst%SpcArr3D, Inst%HcoIDs(N),"
)
foreach(TEXT IN LISTS REQUIRED_GFED_TEXT)
    string(FIND "${GFED_SOURCE}" "${TEXT}" POSITION)
    if(POSITION EQUAL -1)
        message(FATAL_ERROR "Missing GFED scaling integration: ${TEXT}")
    endif()
endforeach()

# A 2-D delivery would leave requested 3-D Ext111 diagnostics unfilled,
# even though the model's actual surface emissions are correct.
string(FIND "${GFED_SOURCE}"
    "CALL HCO_EmisAdd( HcoState, SpcArr, Inst%HcoIDs(N)," SURFACE_2D_POSITION)
if(NOT SURFACE_2D_POSITION EQUAL -1)
    message(FATAL_ERROR "GFED surface emissions must also deliver a 3-D profile")
endif()
# The separate FRAC_OF_PBL%DoUse guard must remain conditional: a surface
# profile deliberately needs no PBL input. Runtime closure tests exercise
# zero-injection allocation/delivery rather than forbidding that valid guard.

string(FIND "${GFED_SOURCE}"
    "Inst%SOAPfrac = Inst%SOAPfrac * Inst%SpcScal(N)" STALE_POSITION)
if(NOT STALE_POSITION EQUAL -1)
    message(FATAL_ERROR "Mutable SOAPfrac species-loop scaling remains")
endif()

string(FIND "${SCALING_SOURCE}"
    "SUBROUTINE CONFIGURE_GFED_CO_RATIOS(" HELPER_POSITION)
if(HELPER_POSITION EQUAL -1)
    message(FATAL_ERROR "GFED CO-ratio scaling helper is missing")
endif()

set(REQUIRED_GFED_INCLUDE_TEXT
    "Inst%GFED4_EMFAC(31,1)=Inst%GFED4_EMFAC(1,1)"
    "Inst%GFED4_EMFAC(31,2)=Inst%GFED4_EMFAC(1,2)"
    "Inst%GFED4_EMFAC(31,3)=Inst%GFED4_EMFAC(1,3)"
    "Inst%GFED4_EMFAC(31,4)=Inst%GFED4_EMFAC(1,4)"
    "Inst%GFED4_EMFAC(31,5)=Inst%GFED4_EMFAC(1,5)"
    "Inst%GFED4_EMFAC(31,6)=Inst%GFED4_EMFAC(1,6)"
)
foreach(TEXT IN LISTS REQUIRED_GFED_INCLUDE_TEXT)
    string(FIND "${GFED_INCLUDE_SOURCE}" "${TEXT}" POSITION)
    if(POSITION EQUAL -1)
        message(FATAL_ERROR "SOAP and CO GFED factors differ: ${TEXT}")
    endif()
endforeach()
