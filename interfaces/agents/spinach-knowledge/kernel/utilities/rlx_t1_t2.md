# kernel/utilities/rlx_t1_t2.m

- Signature: `[R1Op,R2Op]=rlx_t1_t2(spin_system,euler_angles)`

## Purpose

Builds separate longitudinal and transverse relaxation superoperators from the extended T1/T2 rates specified for each spin. It supports isotropic scalar rates and orientation-dependent rates supplied as 3-by-3 tensors or function handles.

## Physical / mathematical content

Each basis state receives the sum of the R1 rates for its longitudinal single-spin components and the sum of the R2 rates for its transverse components. Unit-state components contribute no rate. Thus multi-spin orders relax at the sum of the rates of their constituent non-unit single-spin states.

## Numerical / algorithmic content

For a tensor rate, the orientation vector is obtained from the supplied ZYZ active Euler angles in radians and the rate is evaluated as `ort*rate_tensor*ort'`. A function-handle rate is called with the three Euler angles. The routine evaluates these rates for each isotope, then traverses the basis states in a `parfor` loop to form the per-state sums. The outputs are negative diagonal matrices built with `spdiags`.

## Parameters / inputs

- `spin_system` - Spinach system structure containing the `sphten-liouv` basis and per-isotope `rlx.r1_rates` and `rlx.r2_rates` specifications.
- `euler_angles` - three Euler angles in the ZYZ active convention, in radians, specifying orientation relative to the input orientation. Required if any R1 or R2 rate is a 3-by-3 tensor or function handle; has no effect when rates are scalars.

## Outputs

- `R1Op` - diagonal relaxation superoperator containing longitudinal relaxation terms.
- `R2Op` - diagonal relaxation superoperator containing transverse relaxation terms.

## Implementation structure

The routine converts the basis labels to `L` and `M`, evaluates each isotope's R1 and R2 rate as a scalar, tensor contraction, or function-handle call, and rejects unknown rate specifications or complex resulting rates. For each basis state, non-unit spins with `M=0` contribute their R1 rates and spins with `M~=0` contribute their R2 rates. It then returns `-spdiags(r1_diagonal,0,...)` and `-spdiags(r2_diagonal,0,...)`.

## Reference

[Spin Dynamics Wiki: rlx_t1_t2.m](https://spindynamics.org/wiki/index.php?title=rlx_t1_t2.m)
