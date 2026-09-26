# kernel/summaries/summary_mode_mods.m

- Signature: `summary_mode_mods(spin_system,header)`

## Purpose

Prints the summary of spin Hamiltonian modulation by bosonic mode coordinates: derivatives of spin-spin coupling tensors and of effective local fields with respect to dimensionless mode displacement coordinates. Syntax: summary_mode_mods(spin_system,header)

## Physical / mathematical content

## Numerical / algorithmic content

## Parameters / inputs

- spin_system -Spinach spin system description object
- header -a string of text to precede the summary

## Outputs

- this function prints to the console or to the user-specified
- output via report.m function

## Implementation structure

- Prints the summary of spin Hamiltonian modulation by bosonic
- mode coordinates: derivatives of spin-spin coupling tensors and
- of effective local fields with respect to dimensionless mode
- displacement coordinates. Syntax:
- summary_mode_mods(spin_system,header)
- spin_system -Spinach spin system description object
- header -a string of text to precede the summary
- this function prints to the console or to the user-specified
- output via report.m function
- Check consistency
- Print the summary table
- Print coupling tensor derivatives
