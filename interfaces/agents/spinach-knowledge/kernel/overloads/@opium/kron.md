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

- Kronecker products involving an OPIUM object. Syntax:
- c=kron(a,b)
- a,b -Kronecker operands, can be
- matrices or opia
- c -resulting product
- When both are opia
- Return a bigger opium
- When A is an opium
- Inflate and do the kron
- When B is an opium
- Complain and bomb out
- Never do any enemy a small injury for they are like
