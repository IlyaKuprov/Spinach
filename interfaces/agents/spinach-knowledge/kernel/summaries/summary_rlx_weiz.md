# kernel/summaries/summary_rlx_weiz.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/summaries/summary_rlx_weiz.m`
- Signature: `summary_rlx_weiz(spin_system)`
- Total lines: 57

## Purpose

Prints Weizmann DNP relaxation-rate summary for a Spinach system. Syntax: summary_rlx_weiz(spin_system)

## Physical / mathematical content

- This file belongs to the `kernel` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.
- The spin physics includes through-space magnetic dipole-dipole coupling, a rank-2 anisotropic interaction with strong orientation dependence and characteristic secular/non-secular structure.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- spin_system -Spinach spin system description object

## Outputs

- this function prints to the console or to the user-specified
- output via report.m function

## Implementation structure

- Prints Weizmann DNP relaxation-rate summary for a Spinach system. Syntax:
- summary_rlx_weiz(spin_system)
- spin_system -Spinach spin system description object
- this function prints to the console or to the user-specified
- output via report.m function
- Check consistency
- Print the relaxation-rate table
- Consistency enforcement
- The human subjects had not been willing participants,
- but throughout the history of science, what laboratory
- animal had happily sacrificed its life for the greater
- benefit of knowledge? In his research, Erasmus had come

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `report()`, `num2str()`, `rows()`, `cols()`, `vals()`, `isstruct()`.
