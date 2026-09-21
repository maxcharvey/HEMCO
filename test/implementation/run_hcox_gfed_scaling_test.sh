#!/usr/bin/env bash

set -euo pipefail

this_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
hemco_root=$(git -C "${this_dir}" rev-parse --show-toplevel)
gfed="${hemco_root}/src/Extensions/hcox_gfed_mod.F90"
scaling="${hemco_root}/src/Extensions/hcox_gfed_scaling_mod.F90"
build_dir=${1:?"usage: $0 <configured-HEMCO-build-dir>"}

if [[ ! -f "${build_dir}/CMakeCache.txt" ]]; then
    echo "FAIL: not a configured HEMCO CMake build directory: ${build_dir}" >&2
    exit 2
fi

rg -q --fixed-strings "CALL CONFIGURE_GFED_CO_RATIOS(" "${gfed}"
rg -q --fixed-strings "IF ( TRIM(SpcName) == 'SOAP' ) SpcName = 'CO'" "${gfed}"
rg -q --fixed-strings "SUBROUTINE CONFIGURE_GFED_CO_RATIOS(" "${scaling}"
if rg -q --fixed-strings \
    "Inst%SOAPfrac = Inst%SOAPfrac * Inst%SpcScal(N)" "${gfed}"; then
    echo "FAIL: mutable SOAPfrac species-loop scaling remains" >&2
    exit 1
fi

cmake --build "${build_dir}" --target hcox_gfed_scaling_test
ctest --test-dir "${build_dir}" --output-on-failure \
    -R '^hcox_gfed_(scaling|scaling_integration)$'
