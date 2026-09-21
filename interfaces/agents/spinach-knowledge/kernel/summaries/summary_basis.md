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
