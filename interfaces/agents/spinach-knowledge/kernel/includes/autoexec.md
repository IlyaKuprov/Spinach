# kernel/includes/autoexec.m

- Signature: `(script file)`

## Purpose

This include is executed at the start of create.m, it over- rides all user input. A good use case is forcing polyadic or GPU arithmetic, or some other specific hardware or soft- ware configuration.

## Physical / mathematical content

- Include scripts and shared setup fragments. These files implement tightly scoped runtime setup, parallel profiling, resource guards, or shared kernels included by other Spinach routines.

## Numerical / algorithmic content

## Implementation structure

- This include is executed at the start of create.m, it over-
- rides all user input. A good use case is forcing polyadic
- or GPU arithmetic, or some other specific hardware or soft-
- ware configuration.
- Kill the pointless GPU deprecation warning
- Kill stupid ass figure defaults in R2025a and later
- Do not override user spec
- IK group system settings
- Be careful with GPUs
- 4 workers per GPU are safe
- Do nothing
- This relocates the scratch folder
