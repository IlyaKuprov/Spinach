# etc/textbook/rlx_csa.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/etc/textbook/rlx_csa.m`
- Signature: `[r1,r2]=rlx_csa(B0,isotope,Z,tau_c)`
- Total lines: 72

## Purpose

Redfield theory expressions for CSA relaxation, including contributions from the antisymmetric part. Syntax: [r1,r2]=rlx_csa(B0,isotope,Z,tau_c)

## Physical / mathematical content

- This file belongs to the `etc` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- Chemical-shift anisotropy is present: shielding is treated as a second-rank tensor whose orientation relative to the field or rotor axis modulates line shapes and transfer dynamics.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- B0 -magnet field, Tesla
- isotope -the spins involved, e.g. '15N'
- Z -chemical shift tensor, 3x3
- matrix in ppm
- tau_c -second rank (1/6D) rotational
- correlation time, seconds

## Outputs

- r1 -longitudinal relaxation rate, Hz
- r2 -transverse relaxation rate, Hz
- Note: CSA relaxation rate expressions do not depend on
- the spin quantum number

## Implementation structure

- Redfield theory expressions for CSA relaxation, including
- contributions from the antisymmetric part. Syntax:
- [r1,r2]=rlx_csa(B0,isotope,Z,tau_c)
- B0 -magnet field, Tesla
- isotope -the spins involved, e.g. '15N'
- Z -chemical shift tensor, 3x3
- matrix in ppm
- tau_c -second rank (1/6D) rotational
- correlation time, seconds
- r1 -longitudinal relaxation rate, Hz
- r2 -transverse relaxation rate, Hz
- Note: CSA relaxation rate expressions do not depend on

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `spin()`, `blinv()`, `isscalar()`, `ischar()`.
