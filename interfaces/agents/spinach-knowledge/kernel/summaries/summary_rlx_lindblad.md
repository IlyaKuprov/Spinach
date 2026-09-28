# kernel/summaries/summary_rlx_lindblad.m

- Signature: `summary_rlx_lindblad(spin_system,header)`

## Purpose

Prints per-spin Lindblad relaxation-rate values stored in a Spinach spin system.

## Physical / mathematical content

The table lists each spin's isotope, stored R1 and R2 rates, and label. The function reports these fields; it does not calculate the rates.

## Numerical / algorithmic content

The routine iterates over `spin_system.comp.nspins` and reads `spin_system.rlx.lind_r1_rates` and `lind_r2_rates` for the corresponding spin. It checks that `spin_system` is a structure and `header` is a character string.

## Parameters / inputs

- spin_system - Spinach spin system description object
- header - a string of text to precede the summary

## Outputs

- Prints the rate table through `report.m` to the console or user-specified output.

## Implementation structure

After validation and table headings, one row per spin is emitted with its index, isotope, stored R1 and R2 values, and label.
