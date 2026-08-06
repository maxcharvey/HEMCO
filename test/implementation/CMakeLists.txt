add_executable(hcox_gfed_scaling_test
    hcox_gfed_scaling_test.F90
)
target_link_libraries(hcox_gfed_scaling_test
    PRIVATE HCOX
)

add_test(
    NAME hcox_gfed_scaling
    COMMAND hcox_gfed_scaling_test
)
add_test(
    NAME hcox_gfed_scaling_integration
    COMMAND ${CMAKE_COMMAND}
        -DGFED_FILE=${CMAKE_CURRENT_SOURCE_DIR}/../../src/Extensions/hcox_gfed_mod.F90
        -DGFED_INCLUDE_FILE=${CMAKE_CURRENT_SOURCE_DIR}/../../src/Extensions/hcox_gfed_include_gfed4.H
        -DSCALING_FILE=${CMAKE_CURRENT_SOURCE_DIR}/../../src/Extensions/hcox_gfed_scaling_mod.F90
        -P ${CMAKE_CURRENT_SOURCE_DIR}/check_hcox_gfed_integration.cmake
)

set_tests_properties(
    hcox_gfed_scaling
    hcox_gfed_scaling_integration
    PROPERTIES LABELS "implementation;GFED"
)
