# kernel/spin.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/spin.m`
- Signature: `[gamma,multiplicity]=spin(name)`
- Total lines: 989

## Purpose

Database of multiplicities and magnetogyric ratios for sta- ble and long-lived particles, including spin zero. Syntax: [gamma,multiplicity]=spin(name)

## Physical / mathematical content

- This file belongs to the `kernel` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.
- The effective hardware model is a weakly anharmonic oscillator. Duffing nonlinearity breaks equal level spacing and allows qubit-like addressability within a truncated bosonic ladder.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- name -the name of the isotope, e.g. '15N' or
- '195Pt'; special cases:
- 'G' -ghost spin: gamma=0, mult=1
- 'N', 'M' -neutron, muon
- 'E#' -high-spin electron, # is
- an integer specifying the
- multiplicity
- 'C#' -electromagnetic cavity mo-
- de, # is an integer speci-
- fying the number of popu-
- lation levels
- 'V#' -phonon mode, # is an in-
- teger specifying the num-
- ber of population levels
- 'T#' -transmon, # is an integer
- specifying the number of
- energy levels

## Outputs

- gamma -magnetogyric ratio, rad/(s*Tesla);
- zero for cavities, phonons, and
- transmons
- multiplicity -multiplicity (the number of energy
- or population levels)
- Note: data with no source specified was sourced from Google
- and should be double-checked before running producti-
- on calculations.

## Implementation structure

- Database of multiplicities and magnetogyric ratios for sta-
- ble and long-lived particles, including spin zero. Syntax:
- [gamma,multiplicity]=spin(name)
- name -the name of the isotope, e.g. '15N' or
- '195Pt'; special cases:
- 'G' -ghost spin: gamma=0, mult=1
- 'N', 'M' -neutron, muon
- 'E#' -high-spin electron, # is
- an integer specifying the
- multiplicity
- 'C#' -electromagnetic cavity mo-
- de, # is an integer speci-

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `strcmp()`, `name()`, `regexp()`, `str2double()`, `ischar()`.
