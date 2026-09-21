# experiments/pseudocon/geffect.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/pseudocon/geffect.m`
- Signature: `g=geffect(spin_system,states)`
- Total lines: 125

## Purpose

Effective g-tensor for the user-specified Kramers doublet, computed as described in

## Physical / mathematical content

- Paramagnetic-pseudocontact inference routines. The mathematics includes inverse problems, tensor parameterisation, interpolation, and regularisation.

## Numerical / algorithmic content

- An eigenvalue problem is solved or analysed, so the file is extracting spectra, stationary states, avoided crossings, or modal structure from the effective Hamiltonian or superoperator.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Syntax

```matlab
g=geffect(spin_system,states)
```

## Parameters / inputs

- states -the numbers of the states to use
- (numbered sequentially from the
- lowest to the highest energy)

## Outputs

- g -3x3 g-tensor matrix in Bohr mag-
- neton units

## Implementation structure

- Effective g-tensor for the user-specified Kramers
- doublet, computed as described in
- g=geffect(spin_system,states)
- states -the numbers of the states to use
- (numbered sequentially from the
- lowest to the highest energy)
- g -3x3 g-tensor matrix in Bohr mag-
- neton units
- Check consistency
- Get the g-tensor for each spin
- Get Sx, Sy, Sz operators for each spin
- Get magnetic moment operators

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `gtensorof()`, `operator()`, `hamiltonian()`, `assume()`, `orientation()`, `sqrtm()`, `strcmp()`, `any()`, `states()`.
