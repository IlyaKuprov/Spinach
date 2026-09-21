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
