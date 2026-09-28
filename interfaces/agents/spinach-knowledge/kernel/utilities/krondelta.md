# kernel/utilities/krondelta.m

- Signature: `d=krondelta(a,b)`

## Purpose

Return the Kronecker delta for two integer arguments.

## Physical / mathematical content

The function implements the discrete equality indicator: it is one when the arguments are equal and zero otherwise, returned as a logical value.

## Numerical / algorithmic content

The arguments are checked as real integer scalars. Equality determines the Boolean result.

## Parameters / inputs

- `a` - real integer scalar.
- `b` - real integer scalar.

## Outputs

- `d` - logical scalar, true when `a` equals `b` and false otherwise.

## Implementation structure

After input validation, the function assigns `true()` if `a==b` and `false()` otherwise.
