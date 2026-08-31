# kernel/summaries/summary_basis.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/summaries/summary_basis.m`
- Signature: `summary_basis(spin_system)`
- Total lines: 70

## Purpose

Prints basis-set state summary for a Spinach system. Syntax: summary_basis(spin_system)

## Physical / mathematical content

- This file belongs to the `kernel` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 20-21: Check consistency; implemented by `grumble(spin_system)`.
- Lines 23-24: Get the basis dimension; implemented by `nstates=size(spin_system.bas.basis,1)`.

### Control flow inferred from the code

- Line 25: conditional branch on `nstates > spin_system.tols.basis_hush`.
- Line 31: `for` loop over `n=1:nstates`.
- Line 34: `for` loop over `k=1:spin_system.comp.nspins`.
- Line 40: dispatches on `length(proj)`; cases `1`, `2`.

### Key state/data transformations

- Lines 24: computes `nstates` using `nstates=size(spin_system.bas.basis,1)`.
- Lines 32: computes `current_line` using `current_line=blanks(7+8*spin_system.comp.nspins); spin_number=num2str(n)`.
- Lines 33: computes `current_line(1:length(spin_number))` using `current_line(1:length(spin_number))=spin_number`.
- Lines 35: computes `[L,M]` using `[L,M]=lin2lm(spin_system.bas.basis(n,k))`.
- Lines 36: computes `current_line(7+8*(k-1)+1)` using `current_line(7+8*(k-1)+1)='('`.
- Lines 37: computes `current_line(7+8*(k-1)+2)` using `current_line(7+8*(k-1)+2)=num2str(L)`.
- Lines 38: computes `current_line(7+8*(k-1)+3)` using `current_line(7+8*(k-1)+3)=','`.
- Lines 39: computes `proj` using `proj=num2str(M)`.
- Lines 42: computes `current_line(7+8*(k-1)+4)` using `current_line(7+8*(k-1)+4)=proj`.
- Lines 43: computes `current_line(7+8*(k-1)+5)` using `current_line(7+8*(k-1)+5)=')'`.
- Lines 47: computes `current_line(7+8*(k-1)+6)` using `current_line(7+8*(k-1)+6)=')'`.

### Local helper functions

- Line 61: `grumble()` — `function grumble(spin_system)`. Linux is only free if your time has no value. Jamie Zawinski
  - Representative operation: `if ~isstruct(spin_system)`.
  - Representative operation: `error('spin_system must be a structure.')`.

## Parameters / inputs

- spin_system -Spinach spin system description object

## Outputs

- this function prints to the console or to the user-specified
- output via report.m function

## Implementation structure

- Prints basis-set state summary for a Spinach system. Syntax:
- summary_basis(spin_system)
- spin_system -Spinach spin system description object
- this function prints to the console or to the user-specified
- output via report.m function
- Check consistency
- Get the basis dimension
- Consistency enforcement
- Linux is only free if your time has no value.
- Jamie Zawinski

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `report()`, `num2str()`, `summary()`, `blanks()`, `current_line()`, `lin2lm()`, `proj()`, `isstruct()`.
