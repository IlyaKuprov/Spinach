# etc/textbook/rlx_hfc.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/etc/textbook/rlx_hfc.m`
- Signature: `[r1,r2,rx]=rlx_hfc(B0,HFC,spins,tau_c)`
- Total lines: 139

## Purpose

Redfield theory expressions for hyperfine relaxation and cross- relaxation rates, isotropic tumbling in liquid phase. Syntax: [r1,r2,rx]=rlx_hfc(B0,A,spins,tau_c)

## Physical / mathematical content

- This file belongs to the `etc` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 32-33: Check consistency; implemented by `grumble(B0,HFC,spins,tau_c)`.
- Lines 35-36: Blicharsky invariants and diffusion coefficient; implemented by `[LSqHFC,DSqHFC]=blinv(HFC); r_dif_c=1/(6*tau_c)`.
- Lines 38-39: Multiplicities and spin-squares; implemented by `[~,m_a]=spin(spins{1}); s_a=(m_a-1)/2; Ssq_a=s_a*(s_a+1)`.
- Lines 42-43: Zeeman frequencies; implemented by `omega_a=spin(spins{1})*B0; omega_b=spin(spins{2})*B0`.
- Lines 45-47: Textbook equation, R1 first spin, rank 1 component; implemented by `r1(1)=(1/6)*Ssq_b*LSqHFC*(spden(1,r_dif_c,omega_a)+ spden(1,r_dif_c,omega_a-omega_b))`.
- Lines 49-52: Textbook equation, R1 first spin, rank 2 component; implemented by `r1(1)=r1(1)+(2/27)*Ssq_b*DSqHFC*(3*spden(2,r_dif_c,omega_a)+ 6*spden(2,r_dif_c,omega_a+omega_b)+ 1*spden(2,r_dif_c,omega_a-omega_b))`.
- Lines 54-56: Textbook equation, R1 second spin, rank 1 component; implemented by `r1(2)=(1/6)*Ssq_a*LSqHFC*(spden(1,r_dif_c,omega_b)+ spden(1,r_dif_c,omega_a-omega_b))`.
- Lines 58-61: Textbook equation, R1 second spin, rank 2 component; implemented by `r1(2)=r1(2)+(2/27)*Ssq_a*DSqHFC*(3*spden(2,r_dif_c,omega_b)+ 6*spden(2,r_dif_c,omega_a+omega_b)+ 1*spden(2,r_dif_c,omega_a-omega_b))`.
- Lines 63-66: Textbook equation, R2 first spin, rank 1 component; implemented by `r2(1)=(1/12)*Ssq_b*LSqHFC*(1*spden(1,r_dif_c,omega_a)+ 2*spden(1,r_dif_c,omega_b)+ 1*spden(1,r_dif_c,omega_a-omega_b))`.
- Lines 68-73: Textbook equation, R2 first spin, rank 2 component; implemented by `r2(1)=r2(1)+(1/27)*Ssq_b*DSqHFC*(4*spden(2,r_dif_c,0)+ 3*spden(2,r_dif_c,omega_a)+ 6*spden(2,r_dif_c,omega_b)+ 6*spden(2,r_dif_c,omega_a+omega_b)+ 1*spden(2,r_dif_c,omeg…`.
- Lines 75-78: Textbook equation, R2 second spin, rank 1 component; implemented by `r2(2)=(1/12)*Ssq_a*LSqHFC*(1*spden(1,r_dif_c,omega_b)+ 2*spden(1,r_dif_c,omega_a)+ 1*spden(1,r_dif_c,omega_a-omega_b))`.
- Lines 80-85: Textbook equation, R2 second spin, rank 2 component; implemented by `r2(2)=r2(2)+(1/27)*Ssq_a*DSqHFC*(4*spden(2,r_dif_c,0)+ 3*spden(2,r_dif_c,omega_b)+ 6*spden(2,r_dif_c,omega_a)+ 6*spden(2,r_dif_c,omega_a+omega_b)+ 1*spden(2,r_dif_c,omeg…`.
- Lines 87-88: Textbook equation, longitudinal cross-relaxation rate, rank 1 component; implemented by `rx=-(1/6)*sqrt(Ssq_a*Ssq_b)*LSqHFC*spden(1,r_dif_c,omega_a-omega_b)`.
- Lines 90-92: Textbook equation, longitudinal cross-relaxation rate, rank 2 component; implemented by `rx=rx+(2/27)*sqrt(Ssq_a*Ssq_b)*DSqHFC*(6*spden(2,r_dif_c,omega_a+omega_b)- 1*spden(2,r_dif_c,omega_a-omega_b))`.

### Key state/data transformations

- Lines 36: computes `[LSqHFC,DSqHFC]` using `[LSqHFC,DSqHFC]=blinv(HFC); r_dif_c=1/(6*tau_c)`.
- Lines 39: computes `[~,m_a]` using `[~,m_a]=spin(spins{1}); s_a=(m_a-1)/2; Ssq_a=s_a*(s_a+1)`.
- Lines 40: computes `[~,m_b]` using `[~,m_b]=spin(spins{2}); s_b=(m_b-1)/2; Ssq_b=s_b*(s_b+1)`.
- Lines 43: computes `omega_a` using `omega_a=spin(spins{1})*B0; omega_b=spin(spins{2})*B0`.
- Lines 46-47: computes `r1(1)` using `r1(1)=(1/6)*Ssq_b*LSqHFC*(spden(1,r_dif_c,omega_a)+ spden(1,r_dif_c,omega_a-omega_b))`.
- Lines 55-56: computes `r1(2)` using `r1(2)=(1/6)*Ssq_a*LSqHFC*(spden(1,r_dif_c,omega_b)+ spden(1,r_dif_c,omega_a-omega_b))`.
- Lines 64-66: computes `r2(1)` using `r2(1)=(1/12)*Ssq_b*LSqHFC*(1*spden(1,r_dif_c,omega_a)+ 2*spden(1,r_dif_c,omega_b)+ 1*spden(1,r_dif_c,omega_a-omega_b))`.
- Lines 76-78: computes `r2(2)` using `r2(2)=(1/12)*Ssq_a*LSqHFC*(1*spden(1,r_dif_c,omega_b)+ 2*spden(1,r_dif_c,omega_a)+ 1*spden(1,r_dif_c,omega_a-omega_b))`.
- Lines 88: computes `rx` using `rx=-(1/6)*sqrt(Ssq_a*Ssq_b)*LSqHFC*spden(1,r_dif_c,omega_a-omega_b)`.

### Local helper functions

- Line 97: `grumble()` — `function grumble(B0,HFC,spins,tau_c)`.
  - Representative operation: `if (~isnumeric(B0))||(~isreal(B0))||(~isscalar(B0))`.
  - Representative operation: `error('B0 must be a real number.')`.

## Parameters / inputs

- B0 -magnet field, Tesla
- A -3x3 hyperfine coupling tensor,
- not necessarily symmetric, rad/s
- spins -the spins involved, e.g. {'E','15N'},
- one of those must be an electron
- tau_c -rotational correlation time, seconds

## Outputs

- r1 -two longitudinal relaxation rates, Hz
- r2 -two transverse relaxation rates, Hz
- rx -longitudinal cross-relaxation rate, Hz

## Implementation structure

- Redfield theory expressions for hyperfine relaxation and cross-
- relaxation rates, isotropic tumbling in liquid phase. Syntax:
- [r1,r2,rx]=rlx_hfc(B0,A,spins,tau_c)
- B0 -magnet field, Tesla
- A -3x3 hyperfine coupling tensor,
- not necessarily symmetric, rad/s
- spins -the spins involved, e.g. {'E','15N'},
- one of those must be an electron
- tau_c -rotational correlation time, seconds
- r1 -two longitudinal relaxation rates, Hz
- r2 -two transverse relaxation rates, Hz
- rx -longitudinal cross-relaxation rate, Hz

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `blinv()`, `spin()`, `spden()`, `isscalar()`, `iscell()`, `ischar()`.
