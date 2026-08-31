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

## Code-derived implementation details

### Comment-guided execution stages

- Lines 33-34: Check consistency; implemented by `grumble(B0,isotope,Z,tau_c)`.
- Lines 36-37: Get the carrier frequency; implemented by `omega=B0*(1+trace(1e-6*Z)/3)*spin(isotope)`.
- Lines 39-40: Get Blicharski's invariants; implemented by `[Lsq,Dsq]=blinv(1e-6*Z)`.
- Lines 42-44: Do the textbook math; implemented by `r1=(1/2)*Lsq*omega^2*tau_c/(1+9*tau_c^2*omega^2)+ (2/15)*Dsq*omega^2*tau_c/(1+tau_c^2*omega^2)`.

### Key state/data transformations

- Lines 37: computes `omega` using `omega=B0*(1+trace(1e-6*Z)/3)*spin(isotope)`.
- Lines 40: computes `[Lsq,Dsq]` using `[Lsq,Dsq]=blinv(1e-6*Z)`.
- Lines 43-44: computes `r1` using `r1=(1/2)*Lsq*omega^2*tau_c/(1+9*tau_c^2*omega^2)+ (2/15)*Dsq*omega^2*tau_c/(1+tau_c^2*omega^2)`.
- Lines 45-46: computes `r2` using `r2=(1/4)*Lsq*omega^2*tau_c/(1+9*tau_c^2*omega^2)+ (1/45)*Dsq*omega^2*tau_c*(4+3/(1+tau_c^2*omega^2))`.

### Local helper functions

- Line 51: `grumble()` — `function grumble(B0,isotope,Z,tau_c)`.
  - Representative operation: `if (~isnumeric(B0))||(~isreal(B0))||(~isscalar(B0))`.
  - Representative operation: `error('B0 must be a real number.')`.

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
