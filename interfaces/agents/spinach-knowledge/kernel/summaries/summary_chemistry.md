# kernel/summaries/summary_chemistry.m

- Signature: `summary_chemistry(spin_system)`

## Purpose

Prints a summary of chemical subsystems and exchange in a Spinach system. Syntax: `summary_chemistry(spin_system)`.

## Physical / mathematical content

When the system has multiple chemical subsystems, reports which spins belong to each subsystem and, when present, the inter-subsystem first-order reaction rates. If `chem.flux_rate` is present and nonzero, reports its point-to-point flux rates.

## Numerical / algorithmic content

Lists nonzero off-diagonal entries of the reaction-rate matrix and nonzero entries of the flux-rate matrix.

## Parameters / inputs

- `spin_system` - Spinach spin system description object.

## Outputs

- Prints through `report.m` to the console or the user-specified output.

## Implementation structure

- Checks that the input is a structure, reports subsystem membership when there is more than one chemical subsystem, then reports available reaction rates and flux rates.
