# kernel/utilities/dictum.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/dictum.m`
- Signature: `spin_system=dictum(spin_system,spins,strength)`
- Total lines: 153

## Purpose

Overrides default assumptions about interaction terms surviving rotating frame transformations. Syntax: spin_system=dictum(spin_system,spins,strength)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 38-39: Check consistency; implemented by `grumble(spin_system,spins,strength)`.
- Lines 41-42: Numerical specification, coupling; implemented by `if isnumeric(spins)&&(numel(spins)==2)`.
- Lines 44-47: Report previous assumption; implemented by `report(spin_system,['spins ' num2str(spins(1)) ' (' spin_system.comp.isotopes{spins(1)} '), ' num2str(spins(2)) ' (' spin_system.comp.isotopes{spins(2)} ') coupling assu…`.
- Lines 49-50: Modify the assumption; implemented by `spin_system.inter.coupling.strength{spins(1),spins(2)}=strength`.
- Lines 53-54: Report the new assumption; implemented by `report(spin_system,['changed on user''s request to: ' spin_system.inter.coupling.strength{spins(1),spins(2)}])`.
- Lines 56-57: Numerical specification, Zeeman; implemented by `elseif isnumeric(spins)&&isscalar(spins)`.
- Lines 59-61: Report previous assumption; implemented by `report(spin_system,['spin ' num2str(spins) ' (' spin_system.comp.isotopes{spins} '), Zeeman assumption: ' spin_system.inter.zeeman.strength{spins}])`.
- Lines 63-64: Modify the assumption; implemented by `spin_system.inter.zeeman.strength{spins}=strength`.
- Lines 66-67: Report the new assumption; implemented by `report(spin_system,['changed on user''s request to: ' spin_system.inter.zeeman.strength{spins}])`.
- Lines 69-70: Isotope specification, coupling; implemented by `elseif iscell(spins)&&(numel(spins)==2)`.
- Lines 72-73: Loop over spin pairs; implemented by `for n=find(cellfun(@(x)strcmp(spins{1},x),spin_system.comp.isotopes))`.
- Lines 76-79: Report previous assumption; implemented by `report(spin_system,['spins ' num2str(n) ' (' spin_system.comp.isotopes{n} '), ' num2str(k) ' (' spin_system.comp.isotopes{k} ') coupling assumption: ' spin_system.inter.…`.
- Lines 81-82: Modify the assumption; implemented by `spin_system.inter.coupling.strength{n,k}=strength`.
- Lines 85-86: Report the new assumption; implemented by `report(spin_system,['changed on user''s request to: ' spin_system.inter.coupling.strength{n,k}])`.
- Lines 91-92: Isotope specification, Zeeman; implemented by `elseif iscell(spins)&&isscalar(spins)`.
- Lines 94-95: Loop over spins pairs; implemented by `for n=find(cellfun(@(x)strcmp(spins,x),spin_system.comp.isotopes))`.
- Lines 97-99: Report previous assumption; implemented by `report(spin_system,['spin ' num2str(n) ' (' spin_system.comp.isotopes{n} '), Zeeman assumption: ' spin_system.inter.zeeman.strength{n}])`.
- Lines 101-102: Modify the assumption; implemented by `spin_system.inter.zeeman.strength{n}=strength`.

### Control flow inferred from the code

- Line 42: conditional branch on `isnumeric(spins)&&(numel(spins)==2)`.
- Line 73: `for` loop over `n=find(cellfun(@(x)strcmp(spins{1},x),spin_system.comp.isotopes))`.
- Line 74: `for` loop over `k=find(cellfun(@(x)strcmp(spins{2},x),spin_system.comp.isotopes))`.
- Line 95: `for` loop over `n=find(cellfun(@(x)strcmp(spins,x),spin_system.comp.isotopes))`.

### Key state/data transformations

- Lines 50: computes `spin_system.inter.coupling.strength{spins(1),spins(2)}` using `spin_system.inter.coupling.strength{spins(1),spins(2)}=strength`.
- Lines 51: computes `spin_system.inter.coupling.strength{spins(2),spins(1)}` using `spin_system.inter.coupling.strength{spins(2),spins(1)}=strength`.
- Lines 64: computes `spin_system.inter.zeeman.strength{spins}` using `spin_system.inter.zeeman.strength{spins}=strength`.
- Lines 82: computes `spin_system.inter.coupling.strength{n,k}` using `spin_system.inter.coupling.strength{n,k}=strength`.
- Lines 83: computes `spin_system.inter.coupling.strength{k,n}` using `spin_system.inter.coupling.strength{k,n}=strength`.
- Lines 102: computes `spin_system.inter.zeeman.strength{n}` using `spin_system.inter.zeeman.strength{n}=strength`.

### Local helper functions

- Line 119: `grumble()` — `function grumble(spin_system,spins,strength)`.
  - Representative operation: `if (~isfield(spin_system.inter.coupling,'strength'))|| (~isfield(spin_system.inter.zeeman,'strength'))`.
  - Representative operation: `(~isfield(spin_system.inter.zeeman,'strength'))`.

## Parameters / inputs

- spin_system -Spinach spin system information
- object coming out of assume.m
- spins -a vector with one or two numbers
- or a cell array with one or two
- strings, e.g. [2 4] or {'1H'},
- where one element would cause
- Zeeman interaction assumptions
- to be modified, and two elements
- would cause coupling assumptions
- to be modified.
- strength -new strength specification, see
- the source code of assume.m for
- the available strength specs

## Outputs

- spin_system -updated Spinach spin system in-
- formation object that will be
- used by hamiltonian.m to build
- the Hamiltonian

## Implementation structure

- Overrides default assumptions about interaction terms surviving
- rotating frame transformations. Syntax:
- spin_system=dictum(spin_system,spins,strength)
- spin_system -Spinach spin system information
- object coming out of assume.m
- spins -a vector with one or two numbers
- or a cell array with one or two
- strings, e.g. [2 4] or {'1H'},
- where one element would cause
- Zeeman interaction assumptions
- to be modified, and two elements
- would cause coupling assumptions

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `report()`, `num2str()`, `spins()`, `isscalar()`, `iscell()`, `cellfun()`, `strcmp()`, `isfield()`, `assume()`, `ischar()`, `isvector()`, `any()`, `ismember()`.
