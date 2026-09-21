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
