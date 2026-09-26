# kernel/operators/boson_ortho.m

- Signature: `B=boson_ortho(nlevels)`

## Purpose

Orthogonal bosonic monomials calculated from the bosonic mono- mial basis produced by boson_mono(nlevels). Gram-Schmidt or- thogonalisation is used without normalisation. Syntax: B=boson_ortho(nlevels)

## Physical / mathematical content

- Operator-construction utilities. They build bases and irreducible tensor representations for spin, bosonic, and transition operators.

## Numerical / algorithmic content

## Parameters / inputs

- nlevels -number of bosonic ladder population levels

## Outputs

- B -a cell array of orthogonal bosonic monomials

## Implementation structure

- Orthogonal bosonic monomials calculated from the bosonic mono-
- mial basis produced by boson_mono(nlevels). Gram-Schmidt or-
- thogonalisation is used without normalisation. Syntax:
- B=boson_ortho(nlevels)
- nlevels -number of bosonic ladder population levels
- B -a cell array of orthogonal bosonic monomials
- Check consistency
- Bosonic monomials
- Gram-Schmidt
- Consistency enforcement
- "We're all our own prisons. We are each our own wardens.
- We do our own time. Prison Is In Your Mind."
