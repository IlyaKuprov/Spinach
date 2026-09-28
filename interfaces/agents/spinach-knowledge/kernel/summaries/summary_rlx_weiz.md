# kernel/summaries/summary_rlx_weiz.m

- Signature: `summary_rlx_weiz(spin_system)`

## Purpose

Prints the Weizmann DNP relaxation-rate fields already stored in a Spinach system; it does not calculate the rates.

## Parameters / inputs

- `spin_system` - Spinach spin system structure containing the Weizmann relaxation fields.

## Output

No MATLAB output argument; writes through `report`. The report gives electron and nuclear R1 and R2 rates in Hz, followed by each nonzero inter-nuclear dipolar R1 and R2 entry with its row and column indices, also in Hz.

## Behavior

The function checks that `spin_system` is a structure, prints a Weizmann DNP heading, reads the rates from `spin_system.rlx.weiz_r1e`, `weiz_r2e`, `weiz_r1n`, `weiz_r2n`, `weiz_r1d`, and `weiz_r2d`, then prints a closing separator.
