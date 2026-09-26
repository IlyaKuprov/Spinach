# experiments/pseudocon/geffect.m

- Signature: `g=geffect(spin_system,states)`

## Purpose

Effective g-tensor for the user-specified Kramers doublet, computed as described in

## Physical / mathematical content

- Paramagnetic-pseudocontact inference routines. The mathematics includes inverse problems, tensor parameterisation, interpolation, and regularisation.

## Numerical / algorithmic content

- An eigenvalue problem is solved or analysed, so the file is extracting spectra, stationary states, avoided crossings, or modal structure from the effective Hamiltonian or superoperator.

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
