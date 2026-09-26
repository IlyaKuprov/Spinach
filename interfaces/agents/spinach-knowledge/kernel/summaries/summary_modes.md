# kernel/summaries/summary_modes.m

- Signature: `summary_modes(spin_system,header)`

## Purpose

Prints bosonic mode parameter summary for a Spinach system. Syntax: summary_modes(spin_system,header)

## Physical / mathematical content

- The effective hardware model is a weakly anharmonic oscillator. Duffing nonlinearity breaks equal level spacing and allows qubit-like addressability within a truncated bosonic ladder.

## Numerical / algorithmic content

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
