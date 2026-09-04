# experiments/pseudocon/g2chi.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/pseudocon/g2chi.m`
- Signature: `chi=g2chi(g,T,S)`
- Total lines: 61

## Purpose

Calculates a high-termperature estimate of the magnetic suscep- tibility tensor from the user-supplied g-tensor. Syntax: chi=g2chi(g,T,S)

## Physical / mathematical content

- Paramagnetic-pseudocontact inference routines. The mathematics includes inverse problems, tensor parameterisation, interpolation, and regularisation.

## Numerical / algorithmic content

- The susceptibility tensor is the Curie law prefactor times the g-tensor Gram matrix `g*g.'`, which follows from the Spinach magnetic moment convention `mu=-mu_b*g*S/hbar` and is valid for non-symmetric g-tensors.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 25-26: Check consistency; implemented by `grumble(g,T,S)`.
- Lines 28-29: Fundamental constants; implemented by `mu_b=9.274009994e-24`.
- Lines 33-34: Curie law prefactor; implemented by `prefactor=S*(S+1)*mu_0*(mu_b^2)/(3*k_b*T)`.
- Lines 36-37: Compose the susceptibility tensor from the g-tensor Gram matrix; implemented by `chi=1e30*prefactor*(g*g.')`.

### Key state/data transformations

- Lines 29: computes `mu_b` using `mu_b=9.274009994e-24`.
- Lines 30: computes `mu_0` using `mu_0=4*pi*1e-7`.
- Lines 31: computes `k_b` using `k_b=1.38064852e-23`.
- Lines 34: computes `prefactor` using `prefactor=S*(S+1)*mu_0*(mu_b^2)/(3*k_b*T)`.
- Lines 37: computes `chi` using `chi=1e30*prefactor*(g*g.')`.

### Local helper functions

- Line 42: `grumble()` — `function grumble(g,T,S)`.
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
- Curie law prefactor
- Compose the susceptibility tensor from the g-tensor Gram matrix

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `ismatrix()`, `any()`, `isscalar()`.
