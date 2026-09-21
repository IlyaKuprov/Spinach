# kernel/overloads/@polyadic/full.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/overloads/@polyadic/full.m`
- Signature: `answer=full(p)`
- Total lines: 79

## Purpose

Converts a polyadic representation of a matrix into a full mat- rix. Syntax: answer=full(p) The function opens up all the Kronecker products and uses full arithmetic throughout even if some cores are sparse.

## Physical / mathematical content

- Polyadic tensor-product linear algebra. The emphasis is compressed operator representation, deferred algebra, and efficient Kronecker-structured manipulations.

## Numerical / algorithmic content

## Parameters / inputs

- p -a polyadic object

## Outputs

- answer -a full matrix

## Implementation structure

- Converts a polyadic representation of a matrix into a full mat-
- rix. Syntax:
- answer=full(p)
- The function opens up all the Kronecker products and uses full
- arithmetic throughout even if some cores are sparse.
- p -a polyadic object
- answer -a full matrix
- Process nested polyadics
- Find the core dimensions
- Preallocate the answer
- Loop over the sum
- Compute the polyadic

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `cellfun()`.
