# kernel/summaries/summary_zeeman.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/summaries/summary_zeeman.m`
- Signature: `summary_zeeman(spin_system,header)`
- Total lines: 78

## Purpose

Prints Zeeman interaction tensor summary for a Spinach system. Syntax: summary_zeeman(spin_system,header)

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

- Prints Zeeman interaction tensor summary for a Spinach system. Syntax:
- summary_zeeman(spin_system,header)
- spin_system -Spinach spin system description object
- header -a string of text to precede the summary
- this function prints to the console or to the user-specified
- output via report.m function
- Check consistency
- Print the summary table
- Get the isotropic part
- Get the first and second rank parts
- Do the printing
- Print the break line

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `report()`, `mat2sphten()`, `sphten2mat()`, `pad()`, `num2str()`, `isstruct()`, `ischar()`.
