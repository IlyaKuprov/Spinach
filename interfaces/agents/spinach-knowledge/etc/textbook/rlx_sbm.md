# etc/textbook/rlx_sbm.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/etc/textbook/rlx_sbm.m`
- Signature: `[r1,r2]=rlx_sbm(B0,nucleus,dist,a_iso,e_spin,g_eff,t1e,t2e,tau_r)`
- Total lines: 131

## Purpose

Solomon-Bloembergen-Morgan nuclear relaxation rates due to a paramagnetic centre. Syntax: [r1,r2]=rlx_sbm(B0,nucleus,dist,a_iso,e_spin,g_eff,t1e,t2e,tau_r)

## Physical / mathematical content

- This file belongs to the `etc` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.
- The spin physics includes through-space magnetic dipole-dipole coupling, a rank-2 anisotropic interaction with strong orientation dependence and characteristic secular/non-secular structure.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 40-41: Check consistency; implemented by `grumble(B0,nucleus,dist,a_iso,e_spin,g_eff,t1e,t2e,tau_r)`.
- Lines 43-44: Physical constants; implemented by `mu0=1.25663706127e-6`.
- Lines 48-49: Larmor frequencies; implemented by `omega_i=abs(spin(nucleus)*B0)`.
- Lines 52-53: Effective dipolar correlation times; implemented by `tau_c1=1/(1/tau_r+1/t1e)`.
- Lines 56-57: Dipolar relaxation prefactor; implemented by `spin_sq=e_spin*(e_spin+1)`.
- Lines 61-64: Dipolar longitudinal relaxation rate; implemented by `r1_dd=(2/15)*pref_dd*(3*tau_c1/(1+(omega_i*tau_c1)^2)+ 6*tau_c2/(1+((omega_i+omega_s)*tau_c2)^2)+ tau_c2/(1+((omega_i-omega_s)*tau_c2)^2))`.
- Lines 66-71: Dipolar transverse relaxation rate; implemented by `r2_dd=(1/15)*pref_dd*(4*tau_c1+ 3*tau_c1/(1+(omega_i*tau_c1)^2)+ 6*tau_c2/(1+(omega_s*tau_c2)^2)+ 6*tau_c2/(1+((omega_i+omega_s)*tau_c2)^2)+ tau_c2/(1+((omega_i-omega_s)…`.
- Lines 73-75: Contact longitudinal relaxation rate; implemented by `r1_sc=(2/3)*a_iso^2*spin_sq*t2e/ (1+((omega_i-omega_s)*t2e)^2)`.
- Lines 77-79: Contact transverse relaxation rate; implemented by `r2_sc=(1/3)*a_iso^2*spin_sq*(t1e+t2e/ (1+((omega_i-omega_s)*t2e)^2))`.
- Lines 81-82: Package contributions by mechanism; implemented by `r1=[r1_dd r1_sc]`.

### Key state/data transformations

- Lines 44: computes `mu0` using `mu0=1.25663706127e-6`.
- Lines 45: computes `muB` using `muB=9.2740100657e-24`.
- Lines 46: computes `hbar` using `hbar=6.62607015e-34/(2*pi)`.
- Lines 49: computes `omega_i` using `omega_i=abs(spin(nucleus)*B0)`.
- Lines 50: computes `omega_s` using `omega_s=abs(g_eff*muB*B0/hbar)`.
- Lines 53: computes `tau_c1` using `tau_c1=1/(1/tau_r+1/t1e)`.
- Lines 54: computes `tau_c2` using `tau_c2=1/(1/tau_r+1/t2e)`.
- Lines 57: computes `spin_sq` using `spin_sq=e_spin*(e_spin+1)`.
- Lines 58: computes `dist` using `dist=dist*1e-10`.
- Lines 59: computes `pref_dd` using `pref_dd=(mu0/(4*pi))^2*(spin(nucleus)*g_eff*muB)^2*spin_sq/dist^6`.
- Lines 62-64: computes `r1_dd` using `r1_dd=(2/15)*pref_dd*(3*tau_c1/(1+(omega_i*tau_c1)^2)+ 6*tau_c2/(1+((omega_i+omega_s)*tau_c2)^2)+ tau_c2/(1+((omega_i-omega_s)*tau_c2)^2))`.
- Lines 67-71: computes `r2_dd` using `r2_dd=(1/15)*pref_dd*(4*tau_c1+ 3*tau_c1/(1+(omega_i*tau_c1)^2)+ 6*tau_c2/(1+(omega_s*tau_c2)^2)+ 6*tau_c2/(1+((omega_i+omega_s)*tau_c2)^2)+ tau_c2/(1+((omega_i-omega_s)…`.
- Lines 74-75: computes `r1_sc` using `r1_sc=(2/3)*a_iso^2*spin_sq*t2e/ (1+((omega_i-omega_s)*t2e)^2)`.
- Lines 78-79: computes `r2_sc` using `r2_sc=(1/3)*a_iso^2*spin_sq*(t1e+t2e/ (1+((omega_i-omega_s)*t2e)^2))`.
- Lines 82: computes `r1` using `r1=[r1_dd r1_sc]`.
- Lines 83: computes `r2` using `r2=[r2_dd r2_sc]`.

### Local helper functions

- Line 88: `grumble()` — `function grumble(B0,nucleus,dist,a_iso,e_spin,g_eff,t1e,t2e,tau_r)`.
  - Representative operation: `if (~isnumeric(B0))||(~isreal(B0))||(~isscalar(B0))|| (~isfinite(B0))||(B0<=0)`.
  - Representative operation: `(~isfinite(B0))||(B0<=0)`.

## Parameters / inputs

- B0 -magnet field, Tesla
- nucleus -nuclear isotope, e.g. '1H' or '13C'
- dist -electron-nucleus distance, Angstrom
- a_iso -isotropic hyperfine coupling, rad/s
- e_spin -effective electron spin quantum number
- g_eff -effective electron g-factor
- t1e -longitudinal electron relaxation time, seconds
- t2e -transverse electron relaxation time, seconds
- tau_r -rotational correlation time, seconds

## Outputs

- r1 -longitudinal rates [dipolar contact], Hz
- r2 -transverse rates [dipolar contact], Hz
- The spectral density convention is J(omega,tau)=tau/(1+omega^2*tau^2).

## Implementation structure

- Solomon-Bloembergen-Morgan nuclear relaxation rates due to a
- paramagnetic centre. Syntax:
- [r1,r2]=rlx_sbm(B0,nucleus,dist,a_iso,e_spin,g_eff,t1e,t2e,tau_r)
- B0 -magnet field, Tesla
- nucleus -nuclear isotope, e.g. '1H' or '13C'
- dist -electron-nucleus distance, Angstrom
- a_iso -isotropic hyperfine coupling, rad/s
- e_spin -effective electron spin quantum number
- g_eff -effective electron g-factor
- t1e -longitudinal electron relaxation time, seconds
- t2e -transverse electron relaxation time, seconds
- tau_r -rotational correlation time, seconds

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `spin()`, `isscalar()`, `ischar()`, `regexp()`.
