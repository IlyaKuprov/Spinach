# kernel/summaries/summary_pbc.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/summaries/summary_pbc.m`
- Signature: `summary_pbc(spin_system,header)`
- Total lines: 53

## Purpose

Prints periodic boundary condition vector summary for a Spinach system. Syntax: summary_pbc(spin_system,header)

## Physical / mathematical content

- This file belongs to the `kernel` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- spin_system -Spinach spin system description object
- header -a string of text to precede the summary

## Outputs

- this function prints to the console or to the user-specified
- output via report.m function

## Implementation structure

- Prints periodic boundary condition vector summary for a Spinach system. Syntax:
- summary_pbc(spin_system,header)
- spin_system -Spinach spin system description object
- header -a string of text to precede the summary
- this function prints to the console or to the user-specified
- output via report.m function
- Check consistency
- Print the vector table
- Consistency enforcement
- To anger a conservative, lie to him. To
- anger a liberal, tell him the truth.
- Theodore Roosevelt

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `report()`, `pad()`, `num2str()`, `isstruct()`, `ischar()`.
