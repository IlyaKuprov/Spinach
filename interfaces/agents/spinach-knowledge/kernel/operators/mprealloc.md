# kernel/operators/mprealloc.m

- Signature: `A=mprealloc(spin_system,nnzpc)`

## Purpose

Preallocates an operator in the current basis. Syntax: A=mprealloc(spin_system,nnzpc)

## Physical / mathematical content

- Operator-construction utilities. They build bases and irreducible tensor representations for spin, bosonic, and transition operators.

## Numerical / algorithmic content

## Parameters / inputs

- nnzpc -expected number of non-zeros per column

## Outputs

- A -all-zero sparse matrix of the appropriate
- dimension with room for the specified num-
- ber of non-zeroes

## Implementation structure

- Preallocates an operator in the current basis. Syntax:
- A=mprealloc(spin_system,nnzpc)
- nnzpc - expected number of non-zeros per column
- A - all-zero sparse matrix of the appropriate
- dimension with room for the specified num-
- ber of non-zeroes
- Check consistency
- Do the math
- Create a zero Liouville space matrix operator
- Create a zero Hilbert space matrix operator
- Complain and bomb out
- Consistency enforcement
