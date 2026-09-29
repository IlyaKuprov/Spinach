# kernel/overloads/@polyadic/validate.m

[Source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/overloads/@polyadic/validate.m) · [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=polyadic/validate.m)

- Signature: `validate(p)`

## Purpose

Checks selected structural and dimension invariants of a factorised polyadic; it does not construct or materialise the represented matrix and has no output argument.

## Storage checks

The outer `p.cores` array, `p.prefix`, and `p.suffix` must be cell arrays. Each buffered core term must itself be a cell array, and its entries must pass the numeric check. Prefix and suffix entries are also checked as numeric. For entries that pass those checks and are identified as `polyadic` objects, the implementation calls `validate` recursively.

## Dimension checks

For each core term, the validator multiplies the row counts of its core matrices and separately multiplies their column counts. Every buffered term must have the same resulting row and column totals. If present, the last prefix factor's number of columns must equal the core row total; the first suffix factor's number of rows must equal the core column total. Consecutive factors within each prefix or suffix chain must also have matching inner dimensions. A failed type, cell-structure, or dimension check raises an error.

This verifies factor-chain compatibility without multiplying the factors. It does not return the inferred matrix dimensions.

## Warnings

The function prints a warning message if there are more than 100 buffered core terms, more than 100 prefix factors, or more than 100 suffix factors; these thresholds do not raise errors.
