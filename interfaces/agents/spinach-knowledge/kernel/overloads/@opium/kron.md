# kernel/overloads/@opium/kron.m

- Signature: `c=kron(a,b)`

## Purpose

Kronecker products involving an OPIUM object. Syntax: c=kron(a,b)

## Physical / mathematical content

## Numerical / algorithmic content

## Parameters / inputs

- a,b -Kronecker operands, can be
- matrices or opia

## Outputs

- c -resulting product

## Implementation structure

- Return a larger OPIUM object when both operands are OPIUM objects
- Expand either OPIUM operand to a scaled identity matrix before applying `kron`
- Error if either operand is neither numeric nor an OPIUM object
