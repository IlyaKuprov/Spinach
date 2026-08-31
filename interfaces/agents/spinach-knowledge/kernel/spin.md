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

## Code-derived implementation details

### Comment-guided execution stages

- Lines 53-54: Check consistency; implemented by `grumble(name)`.
- Lines 56-57: Process the spec; implemented by `if strcmp(name(1),'E')&&(~isempty(regexp(name,'^E\d','once')))`.
- Lines 59-60: High-spin electrons are a special case; implemented by `multiplicity=str2double(name(2:end))`.
- Lines 63-64: At least two levels; implemented by `if multiplicity<2`.
- Lines 70-72: Cavity modes, phonon modes, and transmons are special cases, no Zeeman interactions; implemented by `multiplicity=str2double(name(2:end)); gamma=0`.
- Lines 74-75: At least 3 levels; implemented by `if multiplicity<3`.
- Lines 81-82: Other cases come from the database; implemented by `switch name`.

### Control flow inferred from the code

- Line 57: conditional branch on `strcmp(name(1),'E')&&(~isempty(regexp(name,'^E\d','once')))`.
- Line 64: conditional branch on `multiplicity<2`.
- Line 75: conditional branch on `multiplicity<3`.
- Line 82: dispatches on `name`; cases `'G'`, `'E'`, `'N'`, `'M'`, `'1H'`, `'2H'`, `'3H'`, `'3He'`, `'4He'`, `'6Li'`, ….

### Key state/data transformations

- Lines 60: computes `multiplicity` using `multiplicity=str2double(name(2:end))`.
- Lines 61: computes `gamma` using `gamma=-1.76085963023e11`.

### Local helper functions

- Line 970: `grumble()` — `function grumble(name)`. I mean that there is no way to disarm any man except through guilt. Through that which he himself has accepted as guilt. If a man has
  - Representative operation: `if ~ischar(name)`.
  - Representative operation: `error('isotope specification must be a character string.')`.

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
