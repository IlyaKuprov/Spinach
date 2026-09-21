# kernel/summaries/summary_chemistry.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/summaries/summary_chemistry.m`
- Signature: `summary_chemistry(spin_system)`
- Total lines: 78

## Purpose

Prints chemical subsystem and exchange summary for a Spinach system. Syntax: summary_chemistry(spin_system)

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

- Prints chemical subsystem and exchange summary for a Spinach system. Syntax:
- summary_chemistry(spin_system)
- spin_system -Spinach spin system description object
- this function prints to the console or to the user-specified
- output via report.m function
- Check consistency
- Report multiple chemical subsystems
- Report spin system partitioning
- Report first-order reaction rates
- Report flux rates if specified
- Consistency enforcement
- Always code as if the guy who ends up maintaining

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `report()`, `num2str()`, `isfield()`, `Rate()`, `strjust()`, `rows()`, `blanks()`, `cols()`, `vals()`, `isstruct()`.
