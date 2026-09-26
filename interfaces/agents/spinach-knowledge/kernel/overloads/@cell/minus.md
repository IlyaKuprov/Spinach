# kernel/overloads/@cell/minus.m

- Signature: `C=minus(A,B)`

## Purpose

Subtracts cell arrays element-by-element. Syntax: A=minus(A,B)

## Physical / mathematical content

## Numerical / algorithmic content

## Parameters / inputs

- A,B -cell arrays of identical topology

## Outputs

- A -the resulting cell array

## Implementation structure

- Subtracts cell arrays element-by-element. Syntax:
- A=minus(A,B)
- A,B -cell arrays of identical topology
- A -the resulting cell array
- Check consistency
- Decide the topology
- Subtract cell-by-cell
- Subtract from each cell
- Complain and bomb out
- Consistency enforcement
- I can, therefore I am.
- Simone Weil
