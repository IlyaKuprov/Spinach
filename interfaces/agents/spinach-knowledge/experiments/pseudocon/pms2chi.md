# experiments/pseudocon/pms2chi.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/pseudocon/pms2chi.m`
- Signature: `[chi,err]=pms2chi(hfcs,shifts,isotopes)`
- Total lines: 110

## Purpose

Runs a least squares fitting procedure on top of Equation 10 from tensor from DFT hyperfine coupling tensors and experimentally ob- served paramagnetic (contact + pseudocontact) shifts. Syntax: chi=pms2chi(hfcs,shifts,isotopes)

## Physical / mathematical content

- Paramagnetic-pseudocontact inference routines. The mathematics includes inverse problems, tensor parameterisation, interpolation, and regularisation.
- The numerical method is quasi-Newton optimisation: curvature information is approximated from successive step and gradient differences instead of forming exact second derivatives every iteration.
- The optimisation logic is Newton or Newton-like: search directions use first- and second-order local curvature information, usually with regularisation or line-search safeguards.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `lsq_err()`, `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- hfcs -cell array of 3x3 hyperfine coupling tensors
- in Gauss
- shifts -vector of the observed pseudocontact shifts,
- excluding the diamagnetic contribution, ppm
- isotopes -cell array of character strings specifying
- isotopes that exhibit each of the chemical
- shifts supplied, for example {'1H','13C'}

## Outputs

- chi -the fitted magnetic susceptibility tensor,
- cubic Angstroms
- err -least squares error
- Note: Gauss units are used for hyperfine couplings because they do
- not depend on the electron g-tensor.

## Implementation structure

- Runs a least squares fitting procedure on top of Equation 10 from
- tensor from DFT hyperfine coupling tensors and experimentally ob-
- served paramagnetic (contact + pseudocontact) shifts. Syntax:
- chi=pms2chi(hfcs,shifts,isotopes)
- hfcs -cell array of 3x3 hyperfine coupling tensors
- in Gauss
- shifts -vector of the observed pseudocontact shifts,
- excluding the diamagnetic contribution, ppm
- isotopes -cell array of character strings specifying
- isotopes that exhibit each of the chemical
- shifts supplied, for example {'1H','13C'}
- chi -the fitted magnetic susceptibility tensor,

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `optimoptions()`, `fminunc()`, `lsq_err()`, `chi()`, `shifts()`, `hfc2pms()`, `iscell()`, `issymmetric()`, `any()`, `ischar()`.
