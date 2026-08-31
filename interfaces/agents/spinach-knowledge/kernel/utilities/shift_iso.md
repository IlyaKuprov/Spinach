# kernel/utilities/shift_iso.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/shift_iso.m`
- Signature: `tensors=shift_iso(tensors,spin_numbers,new_iso)`
- Total lines: 75

## Purpose

Replaces the isotropic parts of interaction tensors with user- supplied values. This is useful for correcting DFT calculations, where the anisotropy of the various spin interactions is usual- ly satisfactory, but the isotropic part is not. Syntax: tensors=shift_iso(tensors,spin_numbers,new_iso)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 34-35: Check consistency; implemented by `grumble(tensors,spin_numbers,new_iso)`.
- Lines 37-38: Loop over the tensors; implemented by `for n=1:numel(spin_numbers)`.
- Lines 40-41: Isolate the anisotropy; implemented by `[~,rank1,rank2]=mat2sphten(tensors{spin_numbers(n)})`.
- Lines 43-44: Rebuild with the new isotropic part; implemented by `tensors{spin_numbers(n)}=sphten2mat([],rank1,rank2)+new_iso(n)*eye(3)`.

### Control flow inferred from the code

- Line 38: `for` loop over `n=1:numel(spin_numbers)`.

### Key state/data transformations

- Lines 41: computes `[~,rank1,rank2]` using `[~,rank1,rank2]=mat2sphten(tensors{spin_numbers(n)})`.
- Lines 44: computes `tensors{spin_numbers(n)}` using `tensors{spin_numbers(n)}=sphten2mat([],rank1,rank2)+new_iso(n)*eye(3)`.

### Local helper functions

- Line 51: `grumble()` — `function grumble(tensors,spin_numbers,new_iso)`.
  - Representative operation: `if (~iscell(tensors))||any(any(~cellfun(@isreal,tensors)))|| any(any(~cellfun(@(x)all(size(x)==[3 3]|isempty(x)),tensors)))`.
  - Representative operation: `any(any(~cellfun(@(x)all(size(x)==[3 3]|isempty(x)),tensors)))`.

## Parameters / inputs

- tensors -a cell array of interaction tensors
- as 3x3 matrices
- spin_numbers -a vector containing the numbers
- of spins in the tensors array that
- should have the isotropic parts
- replaced
- new_iso -a vector containing the new isotro-
- pic parts in the same order as the
- spin numbers listed in spin_numbers

## Outputs

- tensors -a cell array of interaction tensors
- as 3x3 matrices

## Implementation structure

- Replaces the isotropic parts of interaction tensors with user-
- supplied values. This is useful for correcting DFT calculations,
- where the anisotropy of the various spin interactions is usual-
- ly satisfactory, but the isotropic part is not. Syntax:
- tensors=shift_iso(tensors,spin_numbers,new_iso)
- tensors -a cell array of interaction tensors
- as 3x3 matrices
- spin_numbers -a vector containing the numbers
- of spins in the tensors array that
- should have the isotropic parts
- replaced
- new_iso -a vector containing the new isotro-

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `mat2sphten()`, `spin_numbers()`, `sphten2mat()`, `new_iso()`, `iscell()`, `any()`, `cellfun()`, `all()`.
