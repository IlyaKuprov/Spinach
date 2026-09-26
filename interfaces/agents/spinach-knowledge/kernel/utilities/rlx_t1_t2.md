# kernel/utilities/rlx_t1_t2.m

- Signature: `[R1Op,R2Op]=rlx_t1_t2(spin_system,euler_angles)`

## Purpose

Extended T1/T2 relaxation model returning the relaxation super- operators separately for the longitudinal and the transverse states. Syntax: [R1Op,R2Op]=rlx_t1_t2(spin_system,euler_angles)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

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
