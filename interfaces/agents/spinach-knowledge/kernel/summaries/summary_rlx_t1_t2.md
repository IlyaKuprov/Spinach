# kernel/summaries/summary_rlx_t1_t2.m

- Signature: `summary_rlx_t1_t2(spin_system,header)`

## Purpose

Prints a table of the stored R1 and R2 relaxation rates for each spin.

## Parameters / inputs

- `spin_system` - Spinach spin system structure.
- `header` - character string printed before the table.

## Output

No MATLAB output argument; writes the header and table through `report`.

## Behavior

The function checks that `spin_system` is a structure and `header` is a character string. It prints one row per spin with its index, isotope, R1 rate, R2 rate, and label. Scalar numeric rates are shown in signed scientific notation; nonscalar numeric rates are labelled “anisotropic”, and nonnumeric rates are labelled “orientation”.
