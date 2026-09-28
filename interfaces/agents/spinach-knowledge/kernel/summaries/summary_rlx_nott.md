# kernel/summaries/summary_rlx_nott.m

- Signature: `summary_rlx_nott(spin_system)`

## Purpose

Prints the Nottingham DNP relaxation-rate values stored in a Spinach spin system.

## Physical / mathematical content

The summary reports electron and nuclear R1 and R2 rates from the Nottingham relaxation fields.

## Numerical / algorithmic content

The function reports `spin_system.rlx.nott_r1e`, `nott_r2e`, `nott_r1n`, and `nott_r2n`, labelled electron or nuclear and printed with Hz units. Its consistency check requires `spin_system` to be a structure.

## Parameters / inputs

- spin_system - Spinach spin system description object

## Outputs

- Prints the Nottingham DNP relaxation summary through `report.m` to the console or user-specified output.

## Implementation structure

After validating the input, the routine prints a heading and the four stored rates: electron R1/R2 and nuclear R1/R2, followed by a separator.
