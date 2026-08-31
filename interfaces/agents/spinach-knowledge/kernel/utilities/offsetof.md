# kernel/utilities/offsetof.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/offsetof.m`
- Signature: `offs=offsetof(spin_system,idx)`
- Total lines: 60

## Purpose

Returns the isotropic Zeeman offset of the specified spin from the pure magnetogyric ratio frequency in the current magnet. Syntax: offs=offsetof(spin_system,idx)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 25-26: Check consistency; implemented by `grumble(spin_system,idx)`.
- Lines 28-29: Pull out the Zeeman tensor; implemented by `offs=spin_system.inter.zeeman.matrix{idx}`.
- Lines 31-32: Subtract the magnet frequency term; implemented by `offs=offs-eye(3)*spin_system.inter.basefrqs(idx)`.
- Lines 34-35: Get the isotropic part in Hz; implemented by `offs=trace(offs)/3; offs=-offs/(2*pi)`.

### Key state/data transformations

- Lines 29: computes `offs` using `offs=spin_system.inter.zeeman.matrix{idx}`.

### Local helper functions

- Line 40: `grumble()` — `function grumble(spin_system,idx)`. At the height of the IRAs terrorist campaign on mainland Britain in December 1974, a bomb was lobbed through the
  - Representative operation: `if (~isnumeric(idx))||(~isreal(idx))|| (~isscalar(idx))||(mod(idx,1)~=0)||(idx<1)`.
  - Representative operation: `(~isscalar(idx))||(mod(idx,1)~=0)||(idx<1)`.

## Parameters / inputs

- idx -index of the spin in sys.isotopes
- array, use idxof() to find index
- by the text label

## Outputs

- offs -offset from the pure magnetogyric
- ratio frequency at the current fi-
- eld, Hz

## Implementation structure

- Returns the isotropic Zeeman offset of the specified spin
- from the pure magnetogyric ratio frequency in the current
- magnet. Syntax:
- offs=offsetof(spin_system,idx)
- idx -index of the spin in sys.isotopes
- array, use idxof() to find index
- by the text label
- offs -offset from the pure magnetogyric
- ratio frequency at the current fi-
- eld, Hz
- Check consistency
- Pull out the Zeeman tensor

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `isscalar()`.
