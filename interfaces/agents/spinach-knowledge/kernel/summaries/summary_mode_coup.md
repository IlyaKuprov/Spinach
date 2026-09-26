# kernel/summaries/summary_mode_coup.m

- Signature: `summary_mode_coup(spin_system,header)`

## Purpose

Prints bosonic mode coupling summary for a Spinach system. This covers mode-mode exchange couplings, cross-Kerr couplings, spin- mode exchange couplings, longitudinal and radiation pressure couplings, and dispersive couplings. Syntax: summary_mode_coup(spin_system,header)

## Physical / mathematical content

## Numerical / algorithmic content

## Parameters / inputs

- spin_system -Spinach spin system description object
- header -a string of text to precede the summary

## Outputs

- this function prints to the console or to the user-specified
- output via report.m function

## Implementation structure

- Prints bosonic mode coupling summary for a Spinach system. This
- covers mode-mode exchange couplings, cross-Kerr couplings, spin-
- mode exchange couplings, longitudinal and radiation pressure
- couplings, and dispersive couplings. Syntax:
- summary_mode_coup(spin_system,header)
- spin_system -Spinach spin system description object
- header -a string of text to precede the summary
- this function prints to the console or to the user-specified
- output via report.m function
- Check consistency
- Print the summary table
- Print exchange couplings
