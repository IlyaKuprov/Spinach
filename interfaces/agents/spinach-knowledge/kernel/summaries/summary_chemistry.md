# kernel/summaries/summary_chemistry.m

- Signature: `summary_chemistry(spin_system)`

## Purpose

Prints chemical subsystem and exchange summary for a Spinach system. Syntax: summary_chemistry(spin_system)

## Physical / mathematical content

## Numerical / algorithmic content

## Parameters / inputs

- spin_system -Spinach spin system description object

## Outputs

- this function prints to the console or to the user-specified
- output via report.m function

## Implementation structure

- Prints chemical subsystem and exchange summary for a Spinach system. Syntax:
- summary_chemistry(spin_system)
- spin_system -Spinach spin system description object
- this function prints to the console or to the user-specified
- output via report.m function
- Check consistency
- Report multiple chemical subsystems
- Report spin system partitioning
- Report first-order reaction rates
- Report flux rates if specified
- Consistency enforcement
- Always code as if the guy who ends up maintaining
