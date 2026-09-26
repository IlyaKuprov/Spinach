# kernel/summaries/summary_couplings.m

- Signature: `summary_couplings(spin_system,header)`

## Purpose

Prints spin-spin coupling tensor summary for a Spinach system. Syntax: summary_couplings(spin_system,header)

## Physical / mathematical content

## Numerical / algorithmic content

## Parameters / inputs

- spin_system -Spinach spin system description object
- header -a string of text to precede the summary

## Outputs

- this function prints to the console or to the user-specified
- output via report.m function

## Implementation structure

- Prints spin-spin coupling tensor summary for a Spinach system. Syntax:
- summary_couplings(spin_system,header)
- spin_system -Spinach spin system description object
- header -a string of text to precede the summary
- this function prints to the console or to the user-specified
- output via report.m function
- Check consistency
- Print the summary table
- Detect significant couplings
- Loop over significant couplings
- Get the isotropic part
- Get the first and second rank parts
