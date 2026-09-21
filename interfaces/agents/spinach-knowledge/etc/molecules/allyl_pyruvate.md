# etc/molecules/allyl_pyruvate.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/etc/molecules/allyl_pyruvate.m`
- Signature: `[sys,inter]=allyl_pyruvate(spins)`
- Total lines: 165

## Purpose

Spin system of allyl pyruvate. Isotropic chemical shifts and J-couplings determined by spectral fitting, coordinates and chemical shift anisotropies from DFT. Syntax: [sys,inter]=allyl_pyruvate(spins)

## Physical / mathematical content

- This file belongs to the `etc` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- spins -a cell array containing the isotopes
- to import, e.g. {'1H','13C'}

## Outputs

- sys, inter -Spinach data structures with the
- specification of the spin system
- Note: 13C-13C J-couplings are not provided -this spin system
- is for natural abundance 13C simulations only.

## Implementation structure

- Spin system of allyl pyruvate. Isotropic chemical shifts and
- J-couplings determined by spectral fitting, coordinates and
- chemical shift anisotropies from DFT. Syntax:
- [sys,inter]=allyl_pyruvate(spins)
- spins -a cell array containing the isotopes
- to import, e.g. {'1H','13C'}
- sys, inter -Spinach data structures with the
- specification of the spin system
- Note: 13C-13C J-couplings are not provided -this spin system
- is for natural abundance 13C simulations only.
- Check consistency
- Spin system

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `idxof()`, `remtrace()`, `ismember()`, `iscell()`, `all()`, `cellfun()`.
