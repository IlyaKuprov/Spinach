# kernel/overloads/@rcv/minus.m

- Signature: `A=minus(A,B)`

## Purpose

Subtracts one RCV object from another. Syntax: A=minus(A,B)

## Physical / mathematical content

- RCV sparse-matrix storage utilities. The focus is data structure design for sparse linear algebra and low-overhead composition of large matrices.

## Numerical / algorithmic content

## Parameters / inputs

- A -left operand
- B -right operand

## Outputs

- A -result A-B as an RCV sparse matrix

## Implementation structure

- Subtracts one RCV object from another. Syntax:
- A=minus(A,B)
- A -left operand
- B -right operand
- A -result A-B as an RCV sparse matrix
- Check consistency
- Just call plus
- Consistency enforcement
- "My cat had been suffering from severe illness over the past month
- or so. This had meant that he had needed increasingly hands-on
- care. Due to a terminal diagnosis the decision to put him to
- sleep was made; that took place on Monday 11th April.
