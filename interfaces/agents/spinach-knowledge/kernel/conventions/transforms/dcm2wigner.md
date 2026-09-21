# kernel/conventions/transforms/dcm2wigner.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/conventions/transforms/dcm2wigner.m`
- Signature: `D=dcm2wigner(dcm)`
- Total lines: 82

## Purpose

Converts a directional cosine matrix into second-rank Wigner function matrix. Syntax: D=dcm2wigner(dcm)

## Physical / mathematical content

- Convention and tensor-transform utilities. They convert among tensor parameterisations, coordinate systems, and unit systems; the underlying mathematics is linear algebra on rank-2 tensors and rotation representations.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- dcm -directional cosine matrix

## Outputs

- D -matrix of second rank Wigner D functions. Rows
- and columns are sorted by descending ranks:
- [D( 2,2) ... D( 2,-2)
- ... ... ...
- D(-2,2) ... D(-2,-2)]
- Notes: the resulting Wigner matrix is to be used as v=W*v, where v is
- a column vector of irreducible spherical tensor coefficients in
- the following order: T(2,2), T(2,1), T(2,0), T(2,-1), T(2,-2).

## Implementation structure

- Converts a directional cosine matrix into second-rank Wigner function
- matrix. Syntax:
- D=dcm2wigner(dcm)
- dcm -directional cosine matrix
- D -matrix of second rank Wigner D functions. Rows
- and columns are sorted by descending ranks:
- [D( 2,2) ... D( 2,-2)
- ... ... ...
- D(-2,2) ... D(-2,-2)]
- a column vector of irreducible spherical tensor coefficients in
- the following order: T(2,2), T(2,1), T(2,0), T(2,-1), T(2,-2).
- Check consistency

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `dcm()`, `all()`.
