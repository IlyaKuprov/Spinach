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

## Code-derived implementation details

### Comment-guided execution stages

- Lines 32-33: Check consistency; implemented by `grumble(csa,chi,B,T)`.
- Lines 35-36: Fundamental constants; implemented by `mu_0=4*pi*1e-7`.
- Lines 39-40: Keep the second rank in the CSA; implemented by `[~,~,rank2]=mat2sphten(csa)`.
- Lines 43-44: Keep the second rank in the susceptibility; implemented by `[~,~,rank2]=mat2sphten(chi)`.
- Lines 47-48: Compute the RACS; implemented by `racs=1e-30*(B^2/(15*mu_0*k_b*T))*trace(csa*chi)`.

### Key state/data transformations

- Lines 36: computes `mu_0` using `mu_0=4*pi*1e-7`.
- Lines 37: computes `k_b` using `k_b=1.38064852e-23`.
- Lines 40: computes `[~,~,rank2]` using `[~,~,rank2]=mat2sphten(csa)`.
- Lines 41: computes `csa` using `csa=sphten2mat(0,[0 0 0],rank2)`.
- Lines 45: computes `chi` using `chi=sphten2mat(0,[0 0 0],rank2)`.
- Lines 48: computes `racs` using `racs=1e-30*(B^2/(15*mu_0*k_b*T))*trace(csa*chi)`.

### Local helper functions

- Line 53: `grumble()` — `function grumble(csa,chi,B,T)`.
  - Representative operation: `if (~isnumeric(csa))||(~isreal(csa))|| (~ismatrix(csa))||any(size(csa)~=[3 3])`.
  - Representative operation: `(~ismatrix(csa))||any(size(csa)~=[3 3])`.

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
