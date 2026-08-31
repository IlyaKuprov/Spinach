# etc/textbook/rlx_dip.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/etc/textbook/rlx_dip.m`
- Signature: `[r1,r2,rx]=rlx_dip(B0,spins,dist,tau_c)`
- Total lines: 102

## Purpose

Redfield theory expressions for dipolar relaxation and cross- relaxation rates, isotropic tumbling in liquid phase. Syntax: [r1,r2,rx]=rlx_dip(B0,spins,dist,tau_c)

## Physical / mathematical content

- This file belongs to the `etc` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- The spin physics includes through-space magnetic dipole-dipole coupling, a rank-2 anisotropic interaction with strong orientation dependence and characteristic secular/non-secular structure.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 30-31: Check consistency; implemented by `grumble(B0,spins,dist,tau_c)`.
- Lines 33-34: Blicharsky invariant and rotational diffusion coefficient; implemented by `[~,~,~,~,DD]=xyz2dd([0 0 0],[0 0 dist],spins{1},spins{2})`.
- Lines 37-38: Multiplicities and spin-squares; implemented by `[~,m_a]=spin(spins{1}); s_a=(m_a-1)/2; Ssq_a=s_a*(s_a+1)`.
- Lines 41-42: Zeeman frequencies; implemented by `omega_a=spin(spins{1})*B0; omega_b=spin(spins{2})*B0`.
- Lines 44-47: Textbook equation, R1 first spin; implemented by `r1(1)=(2/27)*Ssq_b*DSqDD*(3*spden(2,r_dif_c,omega_a)+ 6*spden(2,r_dif_c,omega_a+omega_b)+ 1*spden(2,r_dif_c,omega_a-omega_b))`.
- Lines 49-52: Textbook equation, R1 second spin; implemented by `r1(2)=(2/27)*Ssq_a*DSqDD*(3*spden(2,r_dif_c,omega_b)+ 6*spden(2,r_dif_c,omega_a+omega_b)+ 1*spden(2,r_dif_c,omega_a-omega_b))`.
- Lines 54-59: Textbook equation, R2 first spin; implemented by `r2(1)=(1/27)*Ssq_b*DSqDD*(4*spden(2,r_dif_c,0)+ 3*spden(2,r_dif_c,omega_a)+ 6*spden(2,r_dif_c,omega_b)+ 6*spden(2,r_dif_c,omega_a+omega_b)+ 1*spden(2,r_dif_c,omega_a-ome…`.
- Lines 61-66: Textbook equation, R2 second spin; implemented by `r2(2)=(1/27)*Ssq_a*DSqDD*(4*spden(2,r_dif_c,0)+ 3*spden(2,r_dif_c,omega_b)+ 6*spden(2,r_dif_c,omega_a)+ 6*spden(2,r_dif_c,omega_a+omega_b)+ 1*spden(2,r_dif_c,omega_a-ome…`.
- Lines 68-70: Textbook equation, longitudinal cross-relaxation rate; implemented by `rx=(2/27)*sqrt(Ssq_a*Ssq_b)*DSqDD*(6*spden(2,r_dif_c,omega_a+omega_b)- 1*spden(2,r_dif_c,omega_a-omega_b))`.

### Key state/data transformations

- Lines 34: computes `[~,~,~,~,DD]` using `[~,~,~,~,DD]=xyz2dd([0 0 0],[0 0 dist],spins{1},spins{2})`.
- Lines 35: computes `[~,DSqDD]` using `[~,DSqDD]=blinv(DD); r_dif_c=1/(6*tau_c)`.
- Lines 38: computes `[~,m_a]` using `[~,m_a]=spin(spins{1}); s_a=(m_a-1)/2; Ssq_a=s_a*(s_a+1)`.
- Lines 39: computes `[~,m_b]` using `[~,m_b]=spin(spins{2}); s_b=(m_b-1)/2; Ssq_b=s_b*(s_b+1)`.
- Lines 42: computes `omega_a` using `omega_a=spin(spins{1})*B0; omega_b=spin(spins{2})*B0`.
- Lines 45-47: computes `r1(1)` using `r1(1)=(2/27)*Ssq_b*DSqDD*(3*spden(2,r_dif_c,omega_a)+ 6*spden(2,r_dif_c,omega_a+omega_b)+ 1*spden(2,r_dif_c,omega_a-omega_b))`.
- Lines 50-52: computes `r1(2)` using `r1(2)=(2/27)*Ssq_a*DSqDD*(3*spden(2,r_dif_c,omega_b)+ 6*spden(2,r_dif_c,omega_a+omega_b)+ 1*spden(2,r_dif_c,omega_a-omega_b))`.
- Lines 55-59: computes `r2(1)` using `r2(1)=(1/27)*Ssq_b*DSqDD*(4*spden(2,r_dif_c,0)+ 3*spden(2,r_dif_c,omega_a)+ 6*spden(2,r_dif_c,omega_b)+ 6*spden(2,r_dif_c,omega_a+omega_b)+ 1*spden(2,r_dif_c,omega_a-ome…`.
- Lines 62-66: computes `r2(2)` using `r2(2)=(1/27)*Ssq_a*DSqDD*(4*spden(2,r_dif_c,0)+ 3*spden(2,r_dif_c,omega_b)+ 6*spden(2,r_dif_c,omega_a)+ 6*spden(2,r_dif_c,omega_a+omega_b)+ 1*spden(2,r_dif_c,omega_a-ome…`.
- Lines 69-70: computes `rx` using `rx=(2/27)*sqrt(Ssq_a*Ssq_b)*DSqDD*(6*spden(2,r_dif_c,omega_a+omega_b)- 1*spden(2,r_dif_c,omega_a-omega_b))`.

### Local helper functions

- Line 75: `grumble()` — `function grumble(B0,spins,dist,tau_c)`.
  - Representative operation: `if (~isnumeric(B0))||(~isreal(B0))||(~isscalar(B0))`.
  - Representative operation: `error('B0 must be a real number.')`.

## Parameters / inputs

- B0 -magnet field, Tesla
- spins -the spins involved, e.g. {'1H','15N'}
- dist -inter-spin distance, Angstrom
- tau_c -rotational correlation time, seconds

## Outputs

- r1 -two longitudinal relaxation rates, Hz
- r2 -two transverse relaxation rates, Hz
- rx -longitudinal cross-relaxation rate, Hz

## Implementation structure

- Redfield theory expressions for dipolar relaxation and cross-
- relaxation rates, isotropic tumbling in liquid phase. Syntax:
- [r1,r2,rx]=rlx_dip(B0,spins,dist,tau_c)
- B0 -magnet field, Tesla
- spins -the spins involved, e.g. {'1H','15N'}
- dist -inter-spin distance, Angstrom
- tau_c -rotational correlation time, seconds
- r1 -two longitudinal relaxation rates, Hz
- r2 -two transverse relaxation rates, Hz
- rx -longitudinal cross-relaxation rate, Hz
- Check consistency
- Blicharsky invariant and rotational diffusion coefficient

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `xyz2dd()`, `blinv()`, `spin()`, `spden()`, `isscalar()`, `iscell()`, `ischar()`.
