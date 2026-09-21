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
