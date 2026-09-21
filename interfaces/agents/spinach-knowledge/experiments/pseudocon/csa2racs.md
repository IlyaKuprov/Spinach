# experiments/pseudocon/csa2racs.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/pseudocon/csa2racs.m`
- Signature: `racs=csa2racs(csa,chi,B,T)`
- Total lines: 76

## Purpose

Calculates a high-termperature estimate of the residual aniso- tropic chemcial shift from user-supplied CSA tensor and magne- tic susceptibility tensor. Syntax: racs=csa2racs(csa,chi,B,T)

## Physical / mathematical content

- Paramagnetic-pseudocontact inference routines. The mathematics includes inverse problems, tensor parameterisation, interpolation, and regularisation.
- Chemical-shift anisotropy is present: shielding is treated as a second-rank tensor whose orientation relative to the field or rotor axis modulates line shapes and transfer dynamics.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- csa -3x3 chemical shift tensor in ppm
- chi -3x3 magnetic susceptibility tensor
- in cubic Angstroms
- T -absolute temperature in Kelvin
- B -magnetic induction in Tesla

## Outputs

- racs -residual anisotropic chemical
- shift in ppm
- The function implements Equation (2) from the paper by Otting
- and company: http://dx.doi.org/10.1021/ja0564259

## Implementation structure

- Calculates a high-termperature estimate of the residual aniso-
- tropic chemcial shift from user-supplied CSA tensor and magne-
- tic susceptibility tensor. Syntax:
- racs=csa2racs(csa,chi,B,T)
- csa -3x3 chemical shift tensor in ppm
- chi -3x3 magnetic susceptibility tensor
- in cubic Angstroms
- T -absolute temperature in Kelvin
- B -magnetic induction in Tesla
- racs -residual anisotropic chemical
- shift in ppm
- The function implements Equation (2) from the paper by Otting

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `mat2sphten()`, `sphten2mat()`, `ismatrix()`, `any()`, `isscalar()`.
