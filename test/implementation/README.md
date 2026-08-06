# Focused HEMCO implementation tests

`run_hcox_gfed_scaling_test.sh` compiles only the production GFED scaling
helper and its focused test, then exercises:

- all 24 relative orderings of CO, NAP, SOAP, and FSOAP;
- fixed-control absolute emissions and final SOAP/CO and FSOAP/CO ratios;
- NAP isolation, CO scalar and scale-field inheritance, distinct fractions,
  zero emissions, repeated initialization, and fail-fast configuration errors;
- static integration guards for SOAP-to-CO mapping and removal of mutable
  `SOAPfrac` scaling from the GFED species loop.

Run from the HEMCO root:

```text
bash test/implementation/run_hcox_gfed_scaling_test.sh
```

Standalone HEMCO CMake also registers `hcox_gfed_scaling` and
`hcox_gfed_scaling_integration` with CTest, so the existing Linux, macOS, and
Windows workflows execute them.

Status: `NOT_EXECUTED_PENDING_AUTHORIZED_COMPUTE_NODE_BUILD_SEQUENCE`.
This focused test does not replace full HCOX compilation or runtime 2-D/3-D
HEMCO closure checks.
