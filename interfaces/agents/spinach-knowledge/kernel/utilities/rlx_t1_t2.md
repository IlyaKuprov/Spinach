# kernel/utilities/rlx_t1_t2.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/rlx_t1_t2.m`
- Signature: `[R1Op,R2Op]=rlx_t1_t2(spin_system,euler_angles)`
- Total lines: 191

## Purpose

Extended T1/T2 relaxation model returning the relaxation super- operators separately for the longitudinal and the transverse states. Syntax: [R1Op,R2Op]=rlx_t1_t2(spin_system,euler_angles)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- euler_angles -three Euler angles (ZYZ active convention
- in radians) specifying system orientation
- relative to the input orientation; requi-
- red when R1 and/or R2 rates had been spe-
- cified as 3x3 tensor or a function handle,
- this argument has no effect for R1 and R2
- rates specified as scalars.

## Outputs

- R1Op -relaxation superoperator containing
- all longitudinal relaxation terms
- R2Op -relaxation superoperator containing
- all transverse relaxation terms
- Note: multi-spin orders relax at the sum of the rates of
- their constituent single-spin orders.

## Implementation structure

- Extended T1/T2 relaxation model returning the relaxation super-
- operators separately for the longitudinal and the transverse
- states. Syntax:
- [R1Op,R2Op]=rlx_t1_t2(spin_system,euler_angles)
- euler_angles -three Euler angles (ZYZ active convention
- in radians) specifying system orientation
- relative to the input orientation; requi-
- red when R1 and/or R2 rates had been spe-
- cified as 3x3 tensor or a function handle,
- this argument has no effect for R1 and R2
- rates specified as scalars.
- R1Op -relaxation superoperator containing

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `lin2lm()`, `isscalar()`, `r1_rates()`, `exist()`, `euler2dcm()`, `euler_angles()`, `current_r1_rate()`, `r2_rates()`, `current_r2_rate()`, `any()`, `logical()`, `local_r1_rates()`, `local_r2_rates()`, `r1_diagonal()`, `r2_diagonal()`.
