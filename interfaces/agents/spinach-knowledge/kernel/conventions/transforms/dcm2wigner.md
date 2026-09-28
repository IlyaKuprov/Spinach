# kernel/conventions/transforms/dcm2wigner.m

- Signature: `D=dcm2wigner(dcm)`

## Purpose

Converts a directional cosine matrix into second-rank Wigner function matrix. Syntax: D=dcm2wigner(dcm)

## Physical / mathematical content
The input is a 3×3 directional cosine matrix representing a rotation. The output is its rank-2 Wigner D matrix, with rows and columns ordered by magnetic index 2, 1, 0, −1, −2. It acts on a column vector of irreducible spherical tensor coefficients in that same order.

## Numerical / algorithmic content
The function derives two complex coefficients, A and B, from entries of the directional cosine matrix. It checks their amplitudes against dcm(3,3), then tests their phases against other matrix entries, negating A and retesting if necessary. It sets Z=|A|²−|B|² and constructs the 5×5 Wigner matrix from polynomial expressions in A, B, their complex conjugates, and Z. Failed amplitude or phase checks raise errors.

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
The main function calls the local `grumble` function before computing the coefficients and assembling D. `grumble` requires a real numeric 3×3 input. It checks orthogonality and determinant against tolerances: deviations above 1e−6 produce warnings, while deviations above 1e−2 produce errors. The main function applies further amplitude and phase self-consistency checks before returning D.
