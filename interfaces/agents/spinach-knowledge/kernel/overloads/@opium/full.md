# kernel/overloads/@opium/full.m

- Signature: `M=full(M)`

## Purpose

Converts an OPIUM object into the full scaled unit matrix that it represents. Syntax: M=full(M)

## Physical / mathematical content

## Numerical / algorithmic content

## Parameters / inputs

- M -an OPIUM object

## Outputs

- M -a full scaled unit matrix
- of appropriate dimension

## Implementation structure

- Build the matrix M=coeff*eye(dim) represented by the OPIUM object
