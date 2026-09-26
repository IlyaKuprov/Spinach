# kernel/optimcon/drifts.m

- Signature: `[drifts,spc_dim]=drifts(spin_system,context,...`

## Purpose

Returns a cell array of drift Liouvillians suitable for the control.drifts variable in ensemble control optimisations.

## Physical / mathematical content

- Optimal-control core routines. These files implement GRAPE-style objective evaluation, quasi-Newton search, line search, regularisation, distortion models, and waveform parameterisations.

## Numerical / algorithmic content

## Syntax

```matlab
[drifts,spc_dim]=drifts(spin_system,context,parameters)
```

## Parameters / inputs

- context -a function handle to Spinach context
- responsible for handling the ensemble
- parameters -parameters required by the context
- assumptions -assumptions required by the context

## Outputs

- drifts -a cell array of Liouvillians format-
- ted as {{La},{Lb},...}, one per en-
- semble member
- spc_dim -dimension of the classical dynamics
- subspace (e.g. rotor grid in MAS)

## Implementation structure

- Returns a cell array of drift Liouvillians suitable for the
- control.drifts variable in ensemble control optimisations.
- [drifts,spc_dim]=drifts(spin_system,context,parameters)
- context -a function handle to Spinach context
- responsible for handling the ensemble
- parameters -parameters required by the context
- assumptions -assumptions required by the context
- drifts -a cell array of Liouvillians format-
- ted as {{La},{Lb},...}, one per en-
- semble member
- spc_dim -dimension of the classical dynamics
- subspace (e.g. rotor grid in MAS)
