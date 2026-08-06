#!/usr/bin/env bash

set -euo pipefail

this_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
hemco_root=$(git -C "${this_dir}" rev-parse --show-toplevel)
gfed="${hemco_root}/src/Extensions/hcox_gfed_mod.F90"
scaling="${hemco_root}/src/Extensions/hcox_gfed_scaling_mod.F90"
tmp_dir=$(mktemp -d)
trap 'rm -rf "${tmp_dir}"' EXIT

rg -q --fixed-strings "CALL CONFIGURE_GFED_CO_RATIOS(" "${gfed}"
rg -q --fixed-strings "IF ( TRIM(SpcName) == 'SOAP' ) SpcName = 'CO'" "${gfed}"
rg -q --fixed-strings "SUBROUTINE CONFIGURE_GFED_CO_RATIOS(" "${scaling}"
if rg -q --fixed-strings \
    "Inst%SOAPfrac = Inst%SOAPfrac * Inst%SpcScal(N)" "${gfed}"; then
    echo "FAIL: mutable SOAPfrac species-loop scaling remains" >&2
    exit 1
fi

"${FC:-gfortran}" -cpp -ffree-line-length-none \
    "${hemco_root}/src/Shared/Headers/hco_precision_mod.F90" \
    "${scaling}" \
    "${this_dir}/hcox_gfed_scaling_test.F90" \
    -J "${tmp_dir}" \
    -o "${tmp_dir}/hcox_gfed_scaling_test"

"${tmp_dir}/hcox_gfed_scaling_test"
