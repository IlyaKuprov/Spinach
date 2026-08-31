# experiments/pseudocon/hfc2pcs.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/pseudocon/hfc2pcs.m`
- Signature: `[pcs,pcs_tensor]=hfc2pcs(A,chi,isotope)`
- Total lines: 74

## Purpose

Converts hyperfine coupling tensors and susceptibility tensors into pseudocontact shifts (contact component is not included) using Equa- tion 10 from http://dx.doi.org/10.1039/C4CP03106G. Syntax: [pcs,pcs_tensor]=hfc2pcs(A,chi,isotope)

## Physical / mathematical content

- Paramagnetic-pseudocontact inference routines. The mathematics includes inverse problems, tensor parameterisation, interpolation, and regularisation.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 31-32: Check consistency; implemented by `grumble(A,chi,isotope)`.
- Lines 34-35: Keep rank 2 to generate only the pseudocontact part; implemented by `[~,~,rank2]=mat2sphten(A); A=sphten2mat(0,[0 0 0],rank2)`.
- Lines 38-39: Fundamental constants; implemented by `gamma_n=spin(isotope)`.
- Lines 43-44: Collect fundamental constants; implemented by `C=10^4*gamma_n*hbar*mu0/(4*pi*(1e-10)^3)`.
- Lines 46-47: Compute the full paramagnetic shift tensor; implemented by `pcs_tensor=1e6*(1/(4*pi))*A*chi/C`.
- Lines 49-50: Compute the isotropic part; implemented by `pcs=trace(pcs_tensor)/3`.

### Key state/data transformations

- Lines 35: computes `[~,~,rank2]` using `[~,~,rank2]=mat2sphten(A); A=sphten2mat(0,[0 0 0],rank2)`.
- Lines 39: computes `gamma_n` using `gamma_n=spin(isotope)`.
- Lines 40: computes `hbar` using `hbar=1.05457173e-34`.
- Lines 41: computes `mu0` using `mu0=4*pi*1e-7`.
- Lines 44: computes `C` using `C=10^4*gamma_n*hbar*mu0/(4*pi*(1e-10)^3)`.
- Lines 47: computes `pcs_tensor` using `pcs_tensor=1e6*(1/(4*pi))*A*chi/C`.
- Lines 50: computes `pcs` using `pcs=trace(pcs_tensor)/3`.

### Local helper functions

- Line 55: `grumble()` — `function grumble(A,chi,isotope)`.
  - Representative operation: `if (~isnumeric(A))||(~isreal(A))|| (~issymmetric(A))||(any(size(A)~=3))`.
  - Representative operation: `(~issymmetric(A))||(any(size(A)~=3))`.

## Parameters / inputs

- A -hyperfine coupling tensor, Gauss
- chi -magnetic susceptibility tensor, Angstrom^3
- isotope -isotope i.e. '1H'

## Outputs

- pcs -isotropic pseudocontact shift, ppm
- pcs_tensor -pseudocontact shift tensor, ppm
- Note: Gauss units are used for hyperfine couplings because they do
- not depend on the electron g-tensor.

## Implementation structure

- Converts hyperfine coupling tensors and susceptibility tensors into
- pseudocontact shifts (contact component is not included) using Equa-
- tion 10 from http://dx.doi.org/10.1039/C4CP03106G. Syntax:
- [pcs,pcs_tensor]=hfc2pcs(A,chi,isotope)
- A -hyperfine coupling tensor, Gauss
- chi -magnetic susceptibility tensor, Angstrom^3
- isotope -isotope i.e. '1H'
- pcs -isotropic pseudocontact shift, ppm
- pcs_tensor -pseudocontact shift tensor, ppm
- Note: Gauss units are used for hyperfine couplings because they do
- not depend on the electron g-tensor.
- Check consistency

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `mat2sphten()`, `sphten2mat()`, `spin()`, `issymmetric()`, `any()`, `ischar()`.
