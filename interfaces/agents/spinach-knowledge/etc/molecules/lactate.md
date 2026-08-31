# etc/molecules/lactate.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/etc/molecules/lactate.m`
- Signature: `[sys,inter]=lactate(spins)`
- Total lines: 79

## Purpose

Spin system of 13C-labelled lactate with the OH protons assumed to be in rapid exchange with water. Syntax: [sys,inter,bas]=lactate(spins)

## Physical / mathematical content

- This file belongs to the `etc` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 26-27: Check consistency; implemented by `grumble(spins)`.
- Lines 29-30: Isotopes; implemented by `sys.isotopes={'13C','13C','13C','1H','1H','1H','1H'}`.
- Lines 33-34: Chemical shifts (approximate); implemented by `inter.zeeman.scalar={182.0 68.4 20.6 4.0 1.2 1.2 1.2}`.
- Lines 36-37: J-couplings (very approximate); implemented by `inter.coupling.scalar=cell(7,7)`.
- Lines 57-59: Prune the arrays; implemented by `mask=ismember(sys.isotopes,spins)| ismember(sys.labels,spins)`.

### Key state/data transformations

- Lines 30: computes `sys.isotopes` using `sys.isotopes={'13C','13C','13C','1H','1H','1H','1H'}`.
- Lines 31: computes `sys.labels` using `sys.labels={'CO','CA','CB','HA','HB1','HB2','HB3'}`.
- Lines 34: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={182.0 68.4 20.6 4.0 1.2 1.2 1.2}`.
- Lines 37: computes `inter.coupling.scalar` using `inter.coupling.scalar=cell(7,7)`.
- Lines 38: computes `inter.coupling.scalar{idxof(sys,'CO'),idxof(sys,'HA')}` using `inter.coupling.scalar{idxof(sys,'CO'),idxof(sys,'HA')}=3.9`.
- Lines 39: computes `inter.coupling.scalar{idxof(sys,'CO'),idxof(sys,'HB1')}` using `inter.coupling.scalar{idxof(sys,'CO'),idxof(sys,'HB1')}=3.7`.
- Lines 40: computes `inter.coupling.scalar{idxof(sys,'CO'),idxof(sys,'HB2')}` using `inter.coupling.scalar{idxof(sys,'CO'),idxof(sys,'HB2')}=3.7`.
- Lines 41: computes `inter.coupling.scalar{idxof(sys,'CO'),idxof(sys,'HB3')}` using `inter.coupling.scalar{idxof(sys,'CO'),idxof(sys,'HB3')}=3.7`.
- Lines 42: computes `inter.coupling.scalar{idxof(sys,'CO'),idxof(sys,'CA')}` using `inter.coupling.scalar{idxof(sys,'CO'),idxof(sys,'CA')}=40.0`.
- Lines 43: computes `inter.coupling.scalar{idxof(sys,'CO'),idxof(sys,'CB')}` using `inter.coupling.scalar{idxof(sys,'CO'),idxof(sys,'CB')}=3.0`.
- Lines 44: computes `inter.coupling.scalar{idxof(sys,'CA'),idxof(sys,'HA')}` using `inter.coupling.scalar{idxof(sys,'CA'),idxof(sys,'HA')}=145.0`.
- Lines 45: computes `inter.coupling.scalar{idxof(sys,'CA'),idxof(sys,'HB1')}` using `inter.coupling.scalar{idxof(sys,'CA'),idxof(sys,'HB1')}=4.4`.
- Lines 46: computes `inter.coupling.scalar{idxof(sys,'CA'),idxof(sys,'HB2')}` using `inter.coupling.scalar{idxof(sys,'CA'),idxof(sys,'HB2')}=4.4`.
- Lines 47: computes `inter.coupling.scalar{idxof(sys,'CA'),idxof(sys,'HB3')}` using `inter.coupling.scalar{idxof(sys,'CA'),idxof(sys,'HB3')}=4.4`.
- Lines 48: computes `inter.coupling.scalar{idxof(sys,'CA'),idxof(sys,'CB')}` using `inter.coupling.scalar{idxof(sys,'CA'),idxof(sys,'CB')}=40.0`.
- Lines 49: computes `inter.coupling.scalar{idxof(sys,'CB'),idxof(sys,'HB1')}` using `inter.coupling.scalar{idxof(sys,'CB'),idxof(sys,'HB1')}=128.0`.
- Lines 50: computes `inter.coupling.scalar{idxof(sys,'CB'),idxof(sys,'HB2')}` using `inter.coupling.scalar{idxof(sys,'CB'),idxof(sys,'HB2')}=128.0`.
- Lines 51: computes `inter.coupling.scalar{idxof(sys,'CB'),idxof(sys,'HB3')}` using `inter.coupling.scalar{idxof(sys,'CB'),idxof(sys,'HB3')}=128.0`.

### Local helper functions

- Line 68: `grumble()` — `function grumble(spins)`. Why a hundred? If I were wrong, one would have been enough. Albert Einstein's response to Nazis
  - Representative operation: `if (~iscell(spins))||(~all(cellfun(@ischar,spins)))`.
  - Representative operation: `error('spins argument must be a cell array of character strings.')`.

## Parameters / inputs

- spins -a cell array containing the isotopes to
- import, e.g. {'1H','13C'} or a list of
- atom labels, e.g. {'CO','CA','HA'}; the
- lists can be mixed, e.g. {'1H','CA'}

## Outputs

- sys -Spinach spin system description structure
- inter -Spinach interaction description structure

## Implementation structure

- Spin system of 13C-labelled lactate with the OH protons
- assumed to be in rapid exchange with water. Syntax:
- [sys,inter,bas]=lactate(spins)
- spins -a cell array containing the isotopes to
- import, e.g. {'1H','13C'} or a list of
- atom labels, e.g. {'CO','CA','HA'}; the
- lists can be mixed, e.g. {'1H','CA'}
- sys -Spinach spin system description structure
- inter -Spinach interaction description structure
- Check consistency
- Isotopes
- Chemical shifts (approximate)

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `idxof()`, `ismember()`, `iscell()`, `all()`, `cellfun()`.
