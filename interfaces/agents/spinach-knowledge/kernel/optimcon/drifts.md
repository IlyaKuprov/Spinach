# kernel/optimcon/drifts.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/optimcon/drifts.m)

## Purpose and syntax

`[drifts,spc_dim]=drifts(spin_system,context,parameters,assumptions)` returns drift Liouvillians for `control.drifts` in ensemble-control optimisations. There are no optional arguments or defaults.

## Inputs

- `spin_system`: Spinach spin system passed to the context.
- `context`: Required function handle for the ensemble context, for example `@powder` or `@singlerot`.
- `parameters`: Parameters consumed by the context; the function sets `parameters.sum_up=0` before calling it.
- `assumptions`: Character string passed to the context.

The source validates that `context` is a function handle and `assumptions` is a character array. It provides no defaults, and does not validate `spin_system`, `parameters`, or the structure returned by the context.

## Outputs and construction

The context is called with the spin system, an evolution-generator capture function, the modified parameters, and assumptions. For each returned ensemble member, the function forms `H+1i*R+1i*K`; if that member has exactly five entries, it adds `1i*systems{n}{5}` as the hydrodynamics contribution. The output is a cell array shaped as `{{La},{Lb},...}`, one single-Liouvillian cell per ensemble member.

`spc_dim` is calculated as `size(H,1)/spin_system.bas.offsets(end)` after the member loop. It represents the classical-dynamics subspace dimension (for example, the rotor grid in MAS); it is a dimension ratio, not a physical unit. The source does not state physical units for the Liouvillian terms.

<https://spindynamics.org/wiki/index.php?title=drifts.m>
