# interfaces/gaussian/brokensymm.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/interfaces/gaussian/brokensymm.m`
- Signature: `J=brokensymm(props_sing,props_trip)`
- Total lines: 62

## Purpose

Exchange coupling estimation from a pair of DFT logs using Yamaguchi equation. The notation is: H=-2J*(Sa.Sb)

## Physical / mathematical content

- Gaussian interfaces. These parse quantum-chemistry output into spin Hamiltonian ingredients such as hyperfine, shielding, or exchange parameters.
- The relevant state manifold is the singlet/triplet decomposition, where permutation symmetry controls selection rules, relaxation susceptibility, and convertibility to ordinary magnetisation.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 29-30: Check consistency; implemented by `grumble(props_sing,props_trip)`.
- Lines 36-37: Convert from Hartree to Hz; implemented by `J=6.57968974479e15*J`.

### Key state/data transformations

- Lines 33: computes `J` using `J=(props_trip.energy-props_sing.energy)/ (`.

### Local helper functions

- Line 42: `grumble()` — `function grumble(props_sing,props_trip)`.
  - Representative operation: `if ~isfield(props_sing,'energy')`.
  - Representative operation: `error('props_sing.energy field is missing')`.

## Syntax

```matlab
J=brokensymm(props_sing,props_trip)
```

## Parameters / inputs

- props_sing -the output of gparse for the singlet
- state of the biradical
- props_trip -the output of gparse for the triplet
- state of the biradical

## Outputs

- J -an order-of-magnitude (really rough)
- estimate of exchange coupling, Hz

## Implementation structure

- Exchange coupling estimation from a pair of DFT logs using
- Yamaguchi equation. The notation is:
- H=-2J*(Sa.Sb)
- J=brokensymm(props_sing,props_trip)
- props_sing -the output of gparse for the singlet
- state of the biradical
- props_trip -the output of gparse for the triplet
- J -an order-of-magnitude (really rough)
- estimate of exchange coupling, Hz
- Check consistency
- Eq 6 in https://doi.org/10.1063/1.5144696
- Convert from Hartree to Hz

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `isfield()`.
