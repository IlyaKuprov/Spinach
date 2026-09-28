# kernel/overloads/@opium/mtimes.m

- Signature: `c=mtimes(a,b)`

## Purpose

Matrix products involving an OPIUM object. Syntax: c=mtimes(a,b)

## Physical / mathematical content

## Numerical / algorithmic content

## Parameters / inputs

- a,b -opia or numerical arrays

## Outputs

- c -multiplication result

## Implementation structure

- Scale an OPIUM operand when multiplied by a scalar
- Check dimensions and perform the matrix product when an OPIUM operand is multiplied by a numeric matrix
- Apply the OPIUM object's coefficient and dimensions when both operands are OPIUM objects
- Error when operands are neither numeric nor OPIUM objects
