# kernel/overloads/@ttclass/isreal.m

- Signature: `answer=isreal(tt)`

## Purpose

Tests whether a tensor-train object's stored coefficients and core arrays are real-valued.

## Numerical / algorithmic content

The function first checks all entries of `tt.coeff`. If they are real, it checks the core array for every train and returns early when a non-real core is found. A non-`ttclass` input raises an error.

## Parameters / inputs

- tt -tensor train object

## Outputs

- answer -logical true when all coefficients and core elements of the tensor train are real

## Implementation structure

The coefficient check uses `all(isreal(tt.coeff))`; core checks use `isreal` on each stored core. Core traversal is skipped if a coefficient is non-real.
