# experiments/pseudocon/hfc2pms.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/pseudocon/hfc2pms.m`
- Signature: `[pms,pms_tensor]=hfc2pms(A,chi,isotope)`
- Total lines: 72

## Purpose

Converts hyperfine coupling tensors and susceptibility tensors into paramagnetic shifts (contact + pseudocontact component) using Equa- tion 10 from http://dx.doi.org/10.1039/C4CP03106G. Syntax: [pms,pms_tensor]=hfc2pms(A,chi,isotope)

## Physical / mathematical content

- Paramagnetic-pseudocontact inference routines. The mathematics includes inverse problems, tensor parameterisation, interpolation, and regularisation.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 31-32: Check consistency; implemented by `grumble(A,chi,isotope)`.
- Lines 34-35: Fundamental constants; implemented by `gamma_n=spin(isotope)`.
- Lines 39-40: Collect fundamental constants; implemented by `C=10^4*gamma_n*hbar*mu0/(4*pi*(1e-10)^3)`.
- Lines 42-43: Compute the full paramagnetic shift tensor; implemented by `pms_tensor=1e6*(1/(4*pi))*A*chi/C`.
- Lines 45-46: Compute the isotropic part; implemented by `pms=trace(pms_tensor)/3`.

### Key state/data transformations

- Lines 35: computes `gamma_n` using `gamma_n=spin(isotope)`.
- Lines 36: computes `hbar` using `hbar=1.05457173e-34`.
- Lines 37: computes `mu0` using `mu0=4*pi*1e-7`.
- Lines 40: computes `C` using `C=10^4*gamma_n*hbar*mu0/(4*pi*(1e-10)^3)`.
- Lines 43: computes `pms_tensor` using `pms_tensor=1e6*(1/(4*pi))*A*chi/C`.
- Lines 46: computes `pms` using `pms=trace(pms_tensor)/3`.

### Local helper functions

- Line 51: `grumble()` — `function grumble(A,chi,isotope)`.
  - Representative operation: `if (~isnumeric(A))||(~isreal(A))|| (~issymmetric(A))||(any(size(A)~=3))`.
  - Representative operation: `(~issymmetric(A))||(any(size(A)~=3))`.

## Parameters / inputs

- A -hyperfine coupling tensor, Gauss
- chi -magnetic susceptibility tensor, Angstrom^3
- isotope -isotope i.e. '1H'

## Outputs

- pms -isotropic paramagnetic shift, ppm
- pms_tensor -paramagnetic shift tensor, ppm
- Note: Gauss units are used for hyperfine couplings because they do
- not depend on the electron g-tensor.

## Implementation structure

- Converts hyperfine coupling tensors and susceptibility tensors into
- paramagnetic shifts (contact + pseudocontact component) using Equa-
- tion 10 from http://dx.doi.org/10.1039/C4CP03106G. Syntax:
- [pms,pms_tensor]=hfc2pms(A,chi,isotope)
- A -hyperfine coupling tensor, Gauss
- chi -magnetic susceptibility tensor, Angstrom^3
- isotope -isotope i.e. '1H'
- pms -isotropic paramagnetic shift, ppm
- pms_tensor -paramagnetic shift tensor, ppm
- Note: Gauss units are used for hyperfine couplings because they do
- not depend on the electron g-tensor.
- Check consistency

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `spin()`, `issymmetric()`, `any()`, `ischar()`.
