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

## Code-derived implementation details

### Comment-guided execution stages

- Lines 63-64: Check consistency; implemented by `grumble(B0,tau_c,isotopes,deltas,coords)`.
- Lines 66-67: Carrier frequencies; implemented by `omegaA=B0*(1+trace(1e-6*deltas{1})/3)*spin(isotopes{1})`.
- Lines 70-71: Spectral power density function; implemented by `J=@(omega)tau_c/(1+omega^2*tau_c^2)`.
- Lines 73-74: Anisotropic parts of Zeeman tensors; implemented by `ZA=1e-6*B0*(deltas{1}-eye(3)*trace(deltas{1})/3)*spin(isotopes{1})`.
- Lines 77-78: Dipole-dipole coupling tensor; implemented by `[~,~,~,~,DD]=xyz2dd(coords{1},coords{2},isotopes{1},isotopes{2})`.
- Lines 80-81: Blicharski invariants and products; implemented by `[~,DsqZA]=blinv(ZA); [~,DsqZB]=blinv(ZB); [~,DsqDD]=blinv(DD)`.
- Lines 84-85: Longitudinal relaxation rates; implemented by `A.r1.csa=(2/15)*DsqZA*J(omegaA)`.
- Lines 92-93: Transverse relaxation rates; implemented by `A.r2.csa=(1/45)*DsqZA*(4*J(0)+3*J(omegaA))`.
- Lines 102-104: TROSY component relaxation rates; implemented by `A.trosy.dd=(1/180)*DsqDD*(4*J(0)+3*J(omegaA)+J(omegaA-omegaB)+ 3*J(omegaB)+6*J(omegaA+omegaB))`.
- Lines 117-118: Cross-relaxation rate; implemented by `X=(1/90)*DsqDD*(6*J(omegaA+omegaB)-J(omegaA-omegaB))`.

### Key state/data transformations

- Lines 67: computes `omegaA` using `omegaA=B0*(1+trace(1e-6*deltas{1})/3)*spin(isotopes{1})`.
- Lines 68: computes `omegaB` using `omegaB=B0*(1+trace(1e-6*deltas{2})/3)*spin(isotopes{2})`.
- Lines 71: computes `J` using `J=@(omega)tau_c/(1+omega^2*tau_c^2)`.
- Lines 74: computes `ZA` using `ZA=1e-6*B0*(deltas{1}-eye(3)*trace(deltas{1})/3)*spin(isotopes{1})`.
- Lines 75: computes `ZB` using `ZB=1e-6*B0*(deltas{2}-eye(3)*trace(deltas{2})/3)*spin(isotopes{2})`.
- Lines 78: computes `[~,~,~,~,DD]` using `[~,~,~,~,DD]=xyz2dd(coords{1},coords{2},isotopes{1},isotopes{2})`.
- Lines 81: computes `[~,DsqZA]` using `[~,DsqZA]=blinv(ZA); [~,DsqZB]=blinv(ZB); [~,DsqDD]=blinv(DD)`.
- Lines 82: computes `[~,X_DD_ZA]` using `[~,X_DD_ZA]=blprod(DD,ZA); [~,X_DD_ZB]=blprod(DD,ZB)`.
- Lines 85: computes `A.r1.csa` using `A.r1.csa=(2/15)*DsqZA*J(omegaA)`.
- Lines 86: computes `A.r1.dd` using `A.r1.dd=(1/90)*DsqDD*(3*J(omegaA)+J(omegaA-omegaB)+6*J(omegaA+omegaB))`.
- Lines 87: computes `A.r1.total` using `A.r1.total=A.r1.csa+A.r1.dd`.
- Lines 88: computes `B.r1.csa` using `B.r1.csa=(2/15)*DsqZB*J(omegaB)`.
- Lines 89: computes `B.r1.dd` using `B.r1.dd=(1/90)*DsqDD*(3*J(omegaB)+J(omegaB-omegaA)+6*J(omegaB+omegaA))`.
- Lines 90: computes `B.r1.total` using `B.r1.total=B.r1.csa+B.r1.dd`.
- Lines 93: computes `A.r2.csa` using `A.r2.csa=(1/45)*DsqZA*(4*J(0)+3*J(omegaA))`.
- Lines 94-95: computes `A.r2.dd` using `A.r2.dd=(1/180)*DsqDD*(4*J(0)+3*J(omegaA)+J(omegaA-omegaB)+ 6*J(omegaB)+6*J(omegaA+omegaB))`.
- Lines 96: computes `A.r2.total` using `A.r2.total=A.r2.csa+A.r2.dd`.
- Lines 97: computes `B.r2.csa` using `B.r2.csa=(1/45)*DsqZB*(4*J(0)+3*J(omegaB))`.

### Local helper functions

- Line 123: `grumble()` — `function grumble(B0,tau_c,isotopes,deltas,coords)`.
  - Representative operation: `if (~isnumeric(B0))||(~isreal(B0))||(~isscalar(B0))||(~isfinite(B0))`.
  - Representative operation: `error('B0 must be a finite real number.')`.

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
