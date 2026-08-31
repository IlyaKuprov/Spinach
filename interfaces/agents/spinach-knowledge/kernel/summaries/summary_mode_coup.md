# kernel/summaries/summary_mode_coup.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/summaries/summary_mode_coup.m`
- Signature: `summary_mode_coup(spin_system,header)`
- Total lines: 84

## Purpose

Prints bosonic mode coupling summary for a Spinach system. This covers mode-mode exchange couplings, cross-Kerr couplings, spin- mode exchange couplings, longitudinal and radiation pressure couplings, and dispersive couplings. Syntax: summary_mode_coup(spin_system,header)

## Physical / mathematical content

- This file belongs to the `kernel` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 25-26: Check consistency; implemented by `grumble(spin_system,header)`.
- Lines 28-29: Print the summary table; implemented by `report(spin_system,header)`.
- Lines 34-35: Print exchange couplings; implemented by `[rows,cols]=find(~cellfun(@isempty,spin_system.inter.modes.exchange))`.
- Lines 42-43: Print cross-Kerr couplings; implemented by `[rows,cols]=find(~cellfun(@isempty,spin_system.inter.modes.kerr))`.
- Lines 50-51: Print longitudinal and radiation pressure couplings; implemented by `[rows,cols]=find(~cellfun(@isempty,spin_system.inter.modes.longitudinal))`.
- Lines 63-64: Print dispersive couplings; implemented by `[rows,cols]=find(~cellfun(@isempty,spin_system.inter.modes.dispersive))`.

### Control flow inferred from the code

- Line 36: `for` loop over `n=1:numel(rows)`.
- Line 44: `for` loop over `n=1:numel(rows)`.
- Line 52: `for` loop over `n=1:numel(rows)`.
- Line 53: conditional branch on `all(ismember(spin_system.comp.types([rows(n) cols(n)]),{'C','V','T'}))`.
- Line 65: `for` loop over `n=1:numel(rows)`.

### Key state/data transformations

- Lines 35: computes `[rows,cols]` using `[rows,cols]=find(~cellfun(@isempty,spin_system.inter.modes.exchange))`.
- Lines 54: computes `coup_name` using `coup_name='radiation pressure'`.

### Local helper functions

- Line 75: `grumble()` — `function grumble(spin_system,header)`.
  - Representative operation: `if ~isstruct(spin_system)`.
  - Representative operation: `error('spin_system must be a structure.')`.

## Parameters / inputs

- spin_system -Spinach spin system description object
- header -a string of text to precede the summary

## Outputs

- this function prints to the console or to the user-specified
- output via report.m function

## Implementation structure

- Prints bosonic mode coupling summary for a Spinach system. This
- covers mode-mode exchange couplings, cross-Kerr couplings, spin-
- mode exchange couplings, longitudinal and radiation pressure
- couplings, and dispersive couplings. Syntax:
- summary_mode_coup(spin_system,header)
- spin_system -Spinach spin system description object
- header -a string of text to precede the summary
- this function prints to the console or to the user-specified
- output via report.m function
- Check consistency
- Print the summary table
- Print exchange couplings

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `report()`, `cellfun()`, `pad()`, `num2str()`, `rows()`, `cols()`, `all()`, `ismember()`, `isstruct()`, `ischar()`.
