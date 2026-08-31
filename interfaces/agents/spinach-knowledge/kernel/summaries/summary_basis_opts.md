# kernel/summaries/summary_basis_opts.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/summaries/summary_basis_opts.m`
- Signature: `summary_basis_opts(spin_system)`
- Total lines: 78

## Purpose

Prints basis-set option summary for a Spinach system. Syntax: summary_basis_opts(spin_system)

## Physical / mathematical content

- This file belongs to the `kernel` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 20-21: Check consistency; implemented by `grumble(spin_system)`.
- Lines 23-24: Report the formalism; implemented by `switch spin_system.bas.formalism`.
- Lines 37-38: Report the approximation; implemented by `if strcmp(spin_system.bas.formalism,'sphten-liouv')`.

### Control flow inferred from the code

- Line 24: dispatches on `spin_system.bas.formalism`; cases `'zeeman-wavef'`, `'zeeman-hilb'`, `'zeeman-liouv'`, `'sphten-liouv'`.
- Line 38: conditional branch on `strcmp(spin_system.bas.formalism,'sphten-liouv')`.
- Line 39: dispatches on `spin_system.bas.approximation`; cases `'IK-0'`, `'IK-1'`, `'IK-2'`, `'IK-DNP'`, `'none'`.

### Local helper functions

- Line 66: `grumble()` — `function grumble(spin_system)`. According to a trade legend, Uhlenbeck and Goudsmit (students of Ehrenfest when they stumbled upon the concept of spin) presented
  - Representative operation: `if ~isstruct(spin_system)`.
  - Representative operation: `error('spin_system must be a structure.')`.

## Parameters / inputs

- spin_system -Spinach spin system description object

## Outputs

- this function prints to the console or to the user-specified
- output via report.m function

## Implementation structure

- Prints basis-set option summary for a Spinach system. Syntax:
- summary_basis_opts(spin_system)
- spin_system -Spinach spin system description object
- this function prints to the console or to the user-specified
- output via report.m function
- Check consistency
- Report the formalism
- Report the approximation
- Consistency enforcement
- According to a trade legend, Uhlenbeck and Goudsmit (students of
- Ehrenfest when they stumbled upon the concept of spin) presented
- it to Ehrenfest and said, in effect "here's our theory, but don't

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `report()`, `strcmp()`, `int2str()`, `num2str()`, `isstruct()`.
