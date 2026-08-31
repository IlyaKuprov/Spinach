# etc/molecules/strychnine.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/etc/molecules/strychnine.m`
- Signature: `[sys,inter]=strychnine(spins)`
- Total lines: 200

## Purpose

Spin system of strychnine. Isotropic chemical shifts and J-couplings are taken from "200 and more NMR experiments: a practical course" by by Berger and Braun. Coordinates taken as those of the major confor- mer of strychnine proposed in http://dx.doi.org/10.1039/C0CC04114A

## Physical / mathematical content

- This file belongs to the `etc` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.
- Chemical-shift anisotropy is present: shielding is treated as a second-rank tensor whose orientation relative to the field or rotor axis modulates line shapes and transfer dynamics.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 33-34: Check consistency; implemented by `grumble(spins)`.
- Lines 36-37: Shorthands for human-readable coupling designations below; implemented by `H1=1; H2=2; H3=3; H4=4; H8=5; H11a=6; H11b=7; H12=8; H13=9; H14=10`.
- Lines 41-43: Hydrogen atoms; implemented by `H_isotopes={'1H','1H','1H','1H','1H','1H','1H','1H','1H','1H','1H', '1H','1H','1H','1H','1H','1H','1H','1H','1H','1H','1H'}`.
- Lines 46-47: Proton chemical shifts; implemented by `H_zeeman=cell(1,numH)`.
- Lines 60-61: Proton J-couplings; implemented by `H_coupling=cell(numH)`.
- Lines 78-79: Proton coordinates; implemented by `H_coordinates=cell(1,numH)`.
- Lines 103-106: Carbon and nitrogen atoms; implemented by `CN_isotopes={'13C','13C','13C','13C','13C','13C','13C','13C','15N','13C', '13C','13C','13C','13C','13C','13C','13C','13C','15N','13C', '13C','13C','13C'}`.
- Lines 109-110: Carbon and nitrogen chemical shifts; implemented by `CN_zeeman=cell(1,numCN)`.
- Lines 124-125: Carbon and nitrogen coordinates; implemented by `CN_coordinates=cell(1,numCN+2)`.
- Lines 150-151: Combine the arrays; implemented by `sys.isotopes=[H_isotopes, CN_isotopes]`.
- Lines 157-158: Add 13C-1H J-couplings; implemented by `inter.coupling.scalar{H3,3+numH}= 159.6`.
- Lines 181-182: Prune the arrays; implemented by `mask=ismember(sys.isotopes,spins)`.

### Key state/data transformations

- Lines 37: computes `H1` using `H1=1; H2=2; H3=3; H4=4; H8=5; H11a=6; H11b=7; H12=8; H13=9; H14=10`.
- Lines 38: computes `H15a` using `H15a=11; H15b=12; H16=13; H17a=14; H17b=15; H18a=16; H18b=17; H20a=18`.
- Lines 39: computes `H20b` using `H20b=19; H22=20; H23a=21; H23b=22`.
- Lines 42-43: computes `H_isotopes` using `H_isotopes={'1H','1H','1H','1H','1H','1H','1H','1H','1H','1H','1H', '1H','1H','1H','1H','1H','1H','1H','1H','1H','1H','1H'}`.
- Lines 44: computes `numH` using `numH=numel(H_isotopes)`.
- Lines 47: computes `H_zeeman` using `H_zeeman=cell(1,numH)`.
- Lines 48: computes `H_zeeman{H3}` using `H_zeeman{H3}= 7.255; H_zeeman{H22}= 5.915`.
- Lines 49: computes `H_zeeman{H2}` using `H_zeeman{H2}= 7.098; H_zeeman{H1}= 7.167`.
- Lines 50: computes `H_zeeman{H4}` using `H_zeeman{H4}= 8.092; H_zeeman{H12}= 4.288`.
- Lines 51: computes `H_zeeman{H23b}` using `H_zeeman{H23b}=4.066; H_zeeman{H23a}=4.148`.
- Lines 52: computes `H_zeeman{H16}` using `H_zeeman{H16}= 3.963; H_zeeman{H8}= 3.860`.
- Lines 53: computes `H_zeeman{H20b}` using `H_zeeman{H20b}=2.745; H_zeeman{H20a}=3.716`.
- Lines 54: computes `H_zeeman{H18a}` using `H_zeeman{H18a}=3.219; H_zeeman{H18b}=2.878`.
- Lines 55: computes `H_zeeman{H13}` using `H_zeeman{H13}= 1.276; H_zeeman{H17a}=1.890`.
- Lines 56: computes `H_zeeman{H17b}` using `H_zeeman{H17b}=1.890; H_zeeman{H11a}=3.132`.
- Lines 57: computes `H_zeeman{H11b}` using `H_zeeman{H11b}=2.670; H_zeeman{H14}= 3.150`.
- Lines 58: computes `H_zeeman{H15b}` using `H_zeeman{H15b}=1.462; H_zeeman{H15a}=2.360`.
- Lines 61: computes `H_coupling` using `H_coupling=cell(numH)`.

### Local helper functions

- Line 191: `grumble()` — `function grumble(spins)`. Patience, n. -a minor form of despair, disguised as a virtue. Ambrose Bierce
  - Representative operation: `if (~iscell(spins))||(~all(cellfun(@ischar,spins)))`.
  - Representative operation: `error('spins argument must be a cell array of character strings.')`.

## Syntax

```matlab
[sys,inter]=strychnine(spins)
```

## Parameters / inputs

- spins -a cell array containing the isotopes to
- import, e.g. {'1H','13C'}

## Outputs

- sys, inter -Spinach input data structures
- Note: 13C-13C J-couplings are not provided -this file is only
- suitable for natural abundance 13C simulations.
- Note: CSA tensors are not provided -relaxation theory treatments
- on top of this file would not account for CSA relaxation.
- Note: only shifts and coordinates are provided for 15N nuclei.

## Implementation structure

- Spin system of strychnine. Isotropic chemical shifts and J-couplings
- are taken from "200 and more NMR experiments: a practical course" by
- by Berger and Braun. Coordinates taken as those of the major confor-
- mer of strychnine proposed in http://dx.doi.org/10.1039/C0CC04114A
- [sys,inter]=strychnine(spins)
- spins -a cell array containing the isotopes to
- import, e.g. {'1H','13C'}
- sys, inter -Spinach input data structures
- Note: 13C-13C J-couplings are not provided -this file is only
- suitable for natural abundance 13C simulations.
- Note: CSA tensors are not provided -relaxation theory treatments
- on top of this file would not account for CSA relaxation.

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `ismember()`, `iscell()`, `all()`, `cellfun()`.
