# kernel/overloads/@polyadic/inflate.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/overloads/@polyadic/inflate.m`
- Signature: `answer=inflate(p)`
- Total lines: 96

## Purpose

Converts a polyadic representation of a matrix into a sparse mat- rix. Syntax: answer=inflate(p) The function opens up all the Kronecker products while preserving the sparse type if some cores are sparse.

## Physical / mathematical content

- Polyadic tensor-product linear algebra. The emphasis is compressed operator representation, deferred algebra, and efficient Kronecker-structured manipulations.

## Numerical / algorithmic content

## Parameters / inputs

- p -a polyadic object

## Outputs

- answer -a sparse matrix, except if prefixes or suffixes
- are full (in that case, a full matrix)

## Implementation structure

- Converts a polyadic representation of a matrix into a sparse mat-
- rix. Syntax:
- answer=inflate(p)
- The function opens up all the Kronecker products while preserving
- the sparse type if some cores are sparse.
- p -a polyadic object
- answer -a sparse matrix, except if prefixes or suffixes
- are full (in that case, a full matrix)
- Process nested polyadics
- Find the core dimensions
- Get index arrays going
- Loop over the sum

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `cellfun()`, `cell2mat()`, `clear()`.
