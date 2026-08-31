# kernel/utilities/dilute.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/dilute.m`
- Signature: `subsystems=dilute(spin_system,isotope,tuples)`
- Total lines: 93

## Purpose

Splits the spin system into several independent subsystems, each containing only one instance of a user specified isotope or iso- tope tuple that is deemed "dilute". All spin system data is upda- ted accordingly. Basis set information, if found, is destroyed.

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 38-39: Default is singles; implemented by `if ~exist('tuples','var'), tuples=1; end`.
- Lines 41-42: Check consistency; implemented by `grumble(isotope,tuples)`.
- Lines 44-46: Inform the user; implemented by `report(spin_system,['treating ' isotope ' as a dilute, ' 'picking tuples of size ' num2str(tuples)])`.
- Lines 48-49: Find out which spins belong to the dilute species; implemented by `dilute_spins=find(cellfun(@(x)strcmp(x,isotope),spin_system.comp.isotopes))`.
- Lines 51-53: Report the quantity; implemented by `report(spin_system,[num2str(numel(dilute_spins)) ' instances of ' isotope ' found in the system.'])`.
- Lines 55-56: Disallow picking more than there is; implemented by `if numel(dilute_spins)<tuples`.
- Lines 60-61: Get the combinations; implemented by `combos=nchoosek(dilute_spins,tuples)`.
- Lines 64-65: Preallocate the answer; implemented by `subsystems=cell(size(combos,1),1)`.
- Lines 67-68: Create new spin systems; implemented by `for n=1:size(combos,1)`.
- Lines 72-73: Inform the user; implemented by `report(spin_system,[num2str(numel(subsystems)) ' spin subsystems returned.'])`.

### Control flow inferred from the code

- Line 39: conditional branch on `~exist('tuples','var'), tuples=1; end`.
- Line 56: conditional branch on `numel(dilute_spins)<tuples`.
- Line 68: `for` loop over `n=1:size(combos,1)`.

### Key state/data transformations

- Lines 49: computes `dilute_spins` using `dilute_spins=find(cellfun(@(x)strcmp(x,isotope),spin_system.comp.isotopes))`.
- Lines 61: computes `combos` using `combos=nchoosek(dilute_spins,tuples)`.
- Lines 65: computes `subsystems` using `subsystems=cell(size(combos,1),1)`.
- Lines 69: computes `subsystems{n}` using `subsystems{n}=kill_spin(spin_system,setdiff(dilute_spins,combos(n,:)))`.

### Local helper functions

- Line 78: `grumble()` — `function grumble(isotope,tuples)`. "Dear Mother, good news today."
  - Representative operation: `if ~ischar(isotope)`.
  - Representative operation: `error('isotope must be a character string.')`.

## Syntax

```matlab
subsystems=dilute(spin_system,isotope,tuples)
```

## Parameters / inputs

- spin_system -spin system object obtained
- from create.m function
- isotope -isotope specification string,
- for example '13C'
- tuples -1 returns all subsystems with
- a single instance of the isoto-
- pe, 2 returns all subsystems
- where two instances of the iso-
- tope are present, etc.

## Outputs

- subsystems -a cell array of spin system
- objects, with each cell cor-
- responding to one of the new-
- ly formed isotopomers.

## Implementation structure

- Splits the spin system into several independent subsystems, each
- containing only one instance of a user specified isotope or iso-
- tope tuple that is deemed "dilute". All spin system data is upda-
- ted accordingly. Basis set information, if found, is destroyed.
- subsystems=dilute(spin_system,isotope,tuples)
- spin_system -spin system object obtained
- from create.m function
- isotope -isotope specification string,
- for example '13C'
- tuples -1 returns all subsystems with
- a single instance of the isoto-
- pe, 2 returns all subsystems

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `exist()`, `grumble()`, `report()`, `num2str()`, `cellfun()`, `strcmp()`, `nchoosek()`, `kill_spin()`, `setdiff()`, `combos()`, `ischar()`, `isscalar()`.
