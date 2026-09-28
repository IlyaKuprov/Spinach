# kernel/optimcon/drifts.m

- Signature: `[drifts,spc_dim]=drifts(spin_system,context,...`

## Purpose

Returns a cell array of drift Liouvillians for `control.drifts` in ensemble control optimisations.

## Syntax

```matlab
[drifts,spc_dim]=drifts(spin_system,context,parameters,assumptions)
```

## Parameters / inputs

- `spin_system` — Spinach spin system.
- `context` — function handle to the Spinach context responsible for the ensemble, such as `@powder` or `@singlerot`.
- `parameters` — parameters required by the context. The function sets `parameters.sum_up=0` to disable ensemble summation.
- `assumptions` — assumptions required by the context, supplied as a character string.

## Outputs

- `drifts` — cell array of Liouvillians formatted as `{{La},{Lb},...}`, one per ensemble member. Each drift combines `H+1i*R+1i*K` and includes the hydrodynamics term, when present.
- `spc_dim` — dimension of the classical dynamics subspace (e.g. the rotor grid in MAS).

<https://spindynamics.org/wiki/index.php?title=drifts.m>
