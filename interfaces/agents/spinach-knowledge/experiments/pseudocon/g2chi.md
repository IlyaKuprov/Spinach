# experiments/pseudocon/g2chi.m

- Signature: `chi=g2chi(g,T,S)`

## Purpose

Calculates a high-termperature estimate of the magnetic suscep- tibility tensor from the user-supplied g-tensor. Syntax: chi=g2chi(g,T,S)

## Physical / mathematical content

- Paramagnetic-pseudocontact inference routines. The mathematics includes inverse problems, tensor parameterisation, interpolation, and regularisation.

## Numerical / algorithmic content

- The susceptibility tensor is the Curie law prefactor times the g-tensor Gram matrix `g*g.'`, which follows from the Spinach magnetic moment convention `mu=-mu_b*g*S/hbar` and is valid for non-symmetric g-tensors.

## Parameters / inputs

- g -3x3 g-tensor matrix in Bohr magneton units
- T -absolute temperature in Kelvin
- S -electron spin (1/2, 1, 3/2, etc.)

## Outputs

- chi -3x3 magnetic susceptibility tensor in cubic Angstrom

## Implementation structure

- Calculates a high-termperature estimate of the magnetic suscep-
- tibility tensor from the user-supplied g-tensor. Syntax:
- chi=g2chi(g,T,S)
- g -3x3 g-tensor matrix in Bohr magneton units
- T -absolute temperature in Kelvin
- S -electron spin (1/2, 1, 3/2, etc.)
- chi -3x3 magnetic susceptibility tensor in cubic Angstrom
- Check consistency
- Fundamental constants
- Curie law prefactor
- Compose the susceptibility tensor from the g-tensor Gram matrix
