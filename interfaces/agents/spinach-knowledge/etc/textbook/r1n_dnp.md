# etc/textbook/r1n_dnp.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/etc/textbook/r1n_dnp.m`
- Signature: `R1n=r1n_dnp(B0,T,g,T1e,T1n_bulk,r,bet)`
- Total lines: 95

## Purpose

A simple model for nuclear longitudnal relaxation rate dies to the presence of an unpaired electron in cryogenic DNP settings. [literature reference goes here]

## Physical / mathematical content

- This file belongs to the `etc` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 46-47: Check consistency; implemented by `grumble(B0,T,g,T1e,T1n_bulk,r,bet)`.
- Lines 49-50: Fundamental constants; implemented by `mu0=4*pi*1e-7`.
- Lines 54-55: Nuclear longitudinal relaxation rate, Hz; implemented by `sech_sq=sech(g*muB*B0/(2*kB*T))^2`.

### Key state/data transformations

- Lines 50: computes `mu0` using `mu0=4*pi*1e-7`.
- Lines 51: computes `muB` using `muB=9.274010e-24`.
- Lines 52: computes `kB` using `kB =1.380649e-23`.
- Lines 55: computes `sech_sq` using `sech_sq=sech(g*muB*B0/(2*kB*T))^2`.
- Lines 56: computes `geom_dd` using `geom_dd=(1-3*cos(bet)^2)/(r/1e10)^3`.
- Lines 57: computes `R1n` using `R1n=(((mu0/(4*pi))*(g*muB/B0)*geom_dd)^2)*sech_sq/T1e+1/T1n_bulk`.

### Local helper functions

- Line 62: `grumble()` — `function grumble(B0,T,g,T1e,T1n_bulk,r,bet)`.
  - Representative operation: `if (~isnumeric(B0))||(~isreal(B0))||(~isscalar(B0))|| (~isfinite(B0))||(B0==0)`.
  - Representative operation: `(~isfinite(B0))||(B0==0)`.

## Syntax

```matlab
R1n=r1n_dnp(B0,T,g,T1e,T1n_bulk,r,bet)
```

## Parameters / inputs

- B0 -main magnet field, Tesla
- T -absolute temperature, Kelvin
- g -electron g-factor, Bohr
- magneton units
- T1e -electron longitudinal rela-
- xation time, seconds
- T1n_bulk -nuclear longitudinal relaxa-
- tion time far away from the
- electron
- r -electron-nuclear distance,
- Angstrom
- bet -angle between the magnet field
- and the electron-nuclear direc-
- tion, radians

## Outputs

- R1n -nuclear relaxation rate, Hz

## Implementation structure

- A simple model for nuclear longitudnal relaxation
- rate dies to the presence of an unpaired electron
- in cryogenic DNP settings.
- [literature reference goes here]
- R1n=r1n_dnp(B0,T,g,T1e,T1n_bulk,r,bet)
- B0 -main magnet field, Tesla
- T -absolute temperature, Kelvin
- g -electron g-factor, Bohr
- magneton units
- T1e -electron longitudinal rela-
- xation time, seconds
- T1n_bulk -nuclear longitudinal relaxa-

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `sech()`, `isscalar()`.
