# kernel/summaries/summary_mode_coup.m

- Signature: `summary_mode_coup(spin_system,header)`

## Purpose

Prints the bosonic-mode coupling terms stored in a Spinach spin system.

## Physical / mathematical content

The table distinguishes mode-mode exchange, cross-Kerr, longitudinal, and dispersive couplings. Longitudinal entries involving only component types C, V, and T are labelled radiation-pressure couplings; other longitudinal entries retain the longitudinal label. Spin-mode exchange is included in the documented scope, but this routine reports entries from the mode-coupling matrices shown below.

## Numerical / algorithmic content

The function scans the exchange, Kerr, longitudinal, and dispersive coupling arrays for nonempty entries. Each amplitude is divided by `2*pi` and displayed in Hz. It checks that `spin_system` is a structure and `header` is a character string.

## Parameters / inputs

- spin_system - Spinach spin system description object
- header - a string of text to precede the summary

## Outputs

- Prints a coupling table through `report.m` to the console or user-specified output.

## Implementation structure

After input validation and table headings, the routine loops over nonempty entries in each coupling array and prints the component indices, coupling type, and amplitude in Hz.
