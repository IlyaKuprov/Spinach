# kernel/summaries/summary_zeeman.m

- Signature: `summary_zeeman(spin_system,header)`

## Purpose

Prints Zeeman interaction tensor summary for a Spinach system. Syntax: summary_zeeman(spin_system,header)

## Physical / mathematical content

## Numerical / algorithmic content

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
