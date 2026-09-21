# etc/textbook/rlx_dd_csa.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/etc/textbook/rlx_dd_csa.m`
- Signature: `[A,B,X]=rlx_dd_csa(B0,tau_c,isotopes,deltas,coords)`
- Total lines: 162

## Purpose

Redfield theory expressions for some relaxation and cross- relaxation rates in a CSA-DD-CSA system with two spin-1/2 particles. Syntax: [A,B,X]=rlx_dd_csa(B0,tau_c,isotopes,deltas,coords)

## Physical / mathematical content

- This file belongs to the `etc` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- The spin physics includes through-space magnetic dipole-dipole coupling, a rank-2 anisotropic interaction with strong orientation dependence and characteristic secular/non-secular structure.
- Chemical-shift anisotropy is present: shielding is treated as a second-rank tensor whose orientation relative to the field or rotor axis modulates line shapes and transfer dynamics.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- B0 -magnet field, Tesla
- tau_c -rotational correlation time, seconds
- isotopes -the spins involved, e.g. {'13C','19F'}
- deltas -a cell array with two symmetric 3x3
- chemical shift tensors in ppm
- coords -a cell array with two 1x3 Cartesi-
- an coordinate vectors in Angstrom

## Outputs

- (A,B).r(1,2).csa -CSA contribution to R1 and R2
- rates of spins A and B
- (A,B).r(1,2).dd -dipolar contribution to R1 and
- R2 rates of spins A and B
- (A,B).r(1,2).total -total R1 and R2 rates
- (A,B).trosy.dd -dipole contribution to the
- transverse relaxation rate
- of TROSY doublet components
- of spins A and B
- (A,B).trosy.csa -CSA contribution to the
- transverse relaxation rate
- of TROSY doublet components
- of spins A and B
- (A,B).trosy.xc -cross-correlation contribu-
- tion to the transverse rela-
- xation rate of TROSY doublet
- components of spins A and B
- (A,B).trosy.total_bro -transverse relaxation rate of
- the broad TROSY doublet com-
- ponent of spins A and B
- (A,B).trosy.total_nar -transverse relaxation rate of
- the narrow TROSY doublet com-
- ponent of spins A and B
- X -longitudinal cross-relaxation rate between
- spins A and B

## Implementation structure

- Redfield theory expressions for some relaxation and cross-
- relaxation rates in a CSA-DD-CSA system with two spin-1/2
- particles. Syntax:
- [A,B,X]=rlx_dd_csa(B0,tau_c,isotopes,deltas,coords)
- B0 -magnet field, Tesla
- tau_c -rotational correlation time, seconds
- isotopes -the spins involved, e.g. {'13C','19F'}
- deltas -a cell array with two symmetric 3x3
- chemical shift tensors in ppm
- coords -a cell array with two 1x3 Cartesi-
- an coordinate vectors in Angstrom
- (A,B).r(1,2).csa -CSA contribution to R1 and R2

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `spin()`, `xyz2dd()`, `blinv()`, `blprod()`, `isscalar()`, `iscell()`, `any()`, `cellfun()`, `isequal()`, `issymmetric()`, `isrow()`.
