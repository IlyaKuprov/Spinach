# kernel/overloads/@ttclass/plus.m

- Signature: `a=plus(a,b)`

## Purpose

Adds tensor train objects by concatenating their buffered trains; it does not immediately recompress the sum.

## Physical / mathematical content

The operands must both be `ttclass` objects with the same number of cores and matching sizes. Their coefficients, cores and tolerances are concatenated, and entries with zero coefficients are removed. If all coefficients are zero, the result is set to `0*unit_like(a)`.

## Numerical / algorithmic content

No core-wise recompression is performed by this method.

## Parameters / inputs

- a - a tensor train object
- b - a tensor train object

## Outputs

- c - a tensor train object

## Implementation structure

- Validate object types and matching dimensions.
- Concatenate coefficients, cores and tolerances.
- Remove zero-coefficient entries, or return a zero train if none remain.
