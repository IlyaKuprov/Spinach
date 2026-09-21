# kernel/optimcon/drifts.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/optimcon/drifts.m`
- Signature: `[drifts,spc_dim]=drifts(spin_system,context,...`
- Total lines: 81

## Purpose

Returns a cell array of drift Liouvillians suitable for the control.drifts variable in ensemble control optimisations.

## Physical / mathematical content

- Optimal-control core routines. These files implement GRAPE-style objective evaluation, quasi-Newton search, line search, regularisation, distortion models, and waveform parameterisations.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `varargin()`, `context()`, `ischar()`.
