# kernel/overloads/@cell/mtimes.m

- Signature: `C=mtimes(A,B)`

## Purpose

Multiplies all entries of a cell array by a user-specified matrix. Syntax: C=mtimes(A,B)

## Physical / mathematical content

## Numerical / algorithmic content

## Parameters / inputs

- A -a matrix or a cell array thereof
- B -a matrix or a cell array thereof

## Outputs

- C -the resulting cell array
- Note: both arguments cannot be cell arrays.

## Implementation structure

- Multiplies all entries of a cell array by a user-specified
- matrix. Syntax:
- C=mtimes(A,B)
- A -a matrix or a cell array thereof
- B -a matrix or a cell array thereof
- C -the resulting cell array
- Note: both arguments cannot be cell arrays.
- Check consistency
- Decide the topology
- Multiply every cell from the left
- Multiply every cell from the right
- Complain and bomb out
