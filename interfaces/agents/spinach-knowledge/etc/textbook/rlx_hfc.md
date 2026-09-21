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
