# etc/molecules/strychnine.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/etc/molecules/strychnine.m`
- Signature: `[sys,inter]=strychnine(spins)`
- Total lines: 200

## Purpose

Spin system of strychnine. Isotropic chemical shifts and J-couplings are taken from "200 and more NMR experiments: a practical course" by Berger and Braun, except the one-bond C18-H18b coupling, which is taken from http://dx.doi.org/10.1016/j.jmr.2014.02.003. Coordinates taken as those of the major conformer of strychnine proposed in http://dx.doi.org/10.1039/C0CC04114A

## Physical / mathematical content

- This file belongs to the `etc` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.
- Chemical-shift anisotropy is present: shielding is treated as a second-rank tensor whose orientation relative to the field or rotor axis modulates line shapes and transfer dynamics.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

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
- Berger and Braun, except the one-bond C18-H18b coupling, which is
- taken from http://dx.doi.org/10.1016/j.jmr.2014.02.003. Coordinates
- taken as those of the major conformer of strychnine proposed in
- http://dx.doi.org/10.1039/C0CC04114A
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
