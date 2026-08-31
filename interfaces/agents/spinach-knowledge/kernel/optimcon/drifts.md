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

## Code-derived implementation details

### Comment-guided execution stages

- Lines 33-34: Check consistency; implemented by `grumble(context,assumptions)`.
- Lines 36-37: Disable ensemble summation; implemented by `parameters.sum_up=0`.
- Lines 39-40: Capture evolution generators; implemented by `gen_cap=@(varargin)varargin(3:end)`.
- Lines 42-43: Call the user-specified context; implemented by `systems=context(spin_system,gen_cap,parameters,assumptions)`.
- Lines 45-46: Pull drift Liouvillians; implemented by `drifts=cell(1,numel(systems))`.
- Lines 49-50: Get Liouvillian components; implemented by `H=systems{n}{1}; R=systems{n}{2}; K=systems{n}{3}`.
- Lines 52-53: Assign drift Liouvillians; implemented by `drifts{n}={H+1i*R+1i*K}`.
- Lines 55-56: Get hydrodynamics if present; implemented by `if numel(systems{n})==5`.
- Lines 62-63: Return classical subspace dimension; implemented by `spc_dim=size(H,1)/size(spin_system.bas.basis,1)`.

### Control flow inferred from the code

- Line 47: `for` loop over `n=1:numel(systems)`.
- Line 56: conditional branch on `numel(systems{n})==5`.

### Key state/data transformations

- Lines 37: computes `parameters.sum_up` using `parameters.sum_up=0`.
- Lines 40: computes `gen_cap` using `gen_cap=@(varargin)varargin(3:end)`.
- Lines 43: computes `systems` using `systems=context(spin_system,gen_cap,parameters,assumptions)`.
- Lines 46: computes `drifts` using `drifts=cell(1,numel(systems))`.
- Lines 50: computes `H` using `H=systems{n}{1}; R=systems{n}{2}; K=systems{n}{3}`.
- Lines 53: computes `drifts{n}` using `drifts{n}={H+1i*R+1i*K}`.
- Lines 57: computes `drifts{n}{1}` using `drifts{n}{1}=drifts{n}{1}+1i*systems{n}{5}`.
- Lines 63: computes `spc_dim` using `spc_dim=size(H,1)/size(spin_system.bas.basis,1)`.

### Local helper functions

- Line 68: `grumble()` — `function grumble(context,assumptions)`. And so first she was crying for a long time, and then she became evil.
  - Representative operation: `if ~isa(context,'function_handle')`.
  - Representative operation: `error('context must be a function handle.')`.

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
