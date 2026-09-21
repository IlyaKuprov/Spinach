# kernel/summaries/summary_modes.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/summaries/summary_modes.m`
- Signature: `summary_modes(spin_system,header)`
- Total lines: 66

## Purpose

Prints bosonic mode parameter summary for a Spinach system. Syntax: summary_modes(spin_system,header)

## Physical / mathematical content

- This file belongs to the `kernel` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.
- The effective hardware model is a weakly anharmonic oscillator. Duffing nonlinearity breaks equal level spacing and allows qubit-like addressability within a truncated bosonic ladder.

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

- Prints bosonic mode parameter summary for a Spinach system. Syntax:
- summary_modes(spin_system,header)
- spin_system -Spinach spin system description object
- header -a string of text to precede the summary
- this function prints to the console or to the user-specified
- output via report.m function
- Check consistency
- Print the summary table
- Translate the type letter into a word
- Do the printing in Hz
- Consistency enforcement

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `report()`, `ismember()`, `pad()`, `num2str()`, `isstruct()`, `ischar()`.
