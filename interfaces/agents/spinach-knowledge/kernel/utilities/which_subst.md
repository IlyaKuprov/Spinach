# kernel/utilities/which_subst.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/which_subst.m`
- Signature: `subst=which_subst(spin_system,spins)`
- Total lines: 62

## Purpose

Finds out which substance hosts the specified spins; throws an error if there is more than one. Syntax: subst=which_subst(spin_system,spins)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 20-21: Check consistency; implemented by `grumble(spin_system,spins)`.
- Lines 23-25: Find the substances hosting specified spins; implemented by `subst_mask=cellfun(@(x)any(ismember(spins,x)), spin_system.chem.parts)`.
- Lines 27-28: Only one substance is permitted; implemented by `if nnz(subst_mask)>1`.
- Lines 34-35: Get substance number; implemented by `subst=find(subst_mask)`.
- Lines 37-38: Confirm that all spins are in the same substance; implemented by `if ~all(ismember(spins,spin_system.chem.parts{subst}))`.

### Control flow inferred from the code

- Line 28: conditional branch on `nnz(subst_mask)>1`.
- Line 38: conditional branch on `~all(ismember(spins,spin_system.chem.parts{subst}))`.

### Key state/data transformations

- Lines 24-25: computes `subst_mask` using `subst_mask=cellfun(@(x)any(ismember(spins,x)), spin_system.chem.parts)`.
- Lines 35: computes `subst` using `subst=find(subst_mask)`.

### Local helper functions

- Line 45: `grumble()` — `function grumble(spin_system,spins)`.
  - Representative operation: `if (~isnumeric(spins))||(~isreal(spins))|| (~isvector(spins))||any(mod(spins,1)~=0)||any(spins<1)`.
  - Representative operation: `(~isvector(spins))||any(mod(spins,1)~=0)||any(spins<1)`.

## Parameters / inputs

- spins -a list of positive integers

## Outputs

- subst -a positive integer

## Implementation structure

- Finds out which substance hosts the specified spins;
- throws an error if there is more than one. Syntax:
- subst=which_subst(spin_system,spins)
- spins -a list of positive integers
- subst -a positive integer
- Check consistency
- Find the substances hosting specified spins
- Only one substance is permitted
- Get substance number
- Confirm that all spins are in the same substance
- Consistency enforcement
- I swear to you that to think too much is

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `cellfun()`, `any()`, `ismember()`, `nnz()`, `all()`, `isvector()`.
