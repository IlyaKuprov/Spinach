# experiments/pseudocon/g2chi.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/pseudocon/g2chi.m`
- Signature: `chi=g2chi(g,T,S)`
- Total lines: 64

## Purpose

Calculates a high-termperature estimate of the magnetic suscep- tibility tensor from the user-supplied g-tensor. Syntax: chi=g2chi(g,T,S)

## Physical / mathematical content

- Paramagnetic-pseudocontact inference routines. The mathematics includes inverse problems, tensor parameterisation, interpolation, and regularisation.

## Numerical / algorithmic content

- An eigenvalue problem is solved or analysed, so the file is extracting spectra, stationary states, avoided crossings, or modal structure from the effective Hamiltonian or superoperator.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 25-26: Check consistency; implemented by `grumble(g,T,S)`.
- Lines 28-29: Fundamental constants; implemented by `mu_b=9.274009994e-24`.
- Lines 33-34: Diagonalise the g-tensor; implemented by `[V,D]=eig(g)`.
- Lines 36-37: Apply the Curie law to the eigenvalues; implemented by `D=S*(S+1)*mu_0*(mu_b^2)*(D.^2)/(3*k_b*T)`.
- Lines 39-40: Compose the susceptibility tensor; implemented by `chi=1e30*V*D*inv(V)`.

### Key state/data transformations

- Lines 29: computes `mu_b` using `mu_b=9.274009994e-24`.
- Lines 30: computes `mu_0` using `mu_0=4*pi*1e-7`.
- Lines 31: computes `k_b` using `k_b=1.38064852e-23`.
- Lines 34: computes `[V,D]` using `[V,D]=eig(g)`.
- Lines 37: computes `D` using `D=S*(S+1)*mu_0*(mu_b^2)*(D.^2)/(3*k_b*T)`.
- Lines 40: computes `chi` using `chi=1e30*V*D*inv(V)`.

### Local helper functions

- Line 45: `grumble()` — `function grumble(g,T,S)`.
  - Representative operation: `if (~isnumeric(g))||(~isreal(g))|| (~ismatrix(g))||any(size(g)~=[3 3])`.
  - Representative operation: `(~ismatrix(g))||any(size(g)~=[3 3])`.

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
- Diagonalise the g-tensor
- Apply the Curie law to the eigenvalues
- Compose the susceptibility tensor

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `inv()`, `ismatrix()`, `any()`, `isscalar()`.
