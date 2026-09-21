# examples/fundamentals/operator_tests/commutation_7.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/fundamentals/operator_tests/commutation_7.m`
- Signature: `commutation_7()`
- Total lines: 126

## Purpose

Expansion relations for operator basis transforms.

## Physical / mathematical content

- Fundamentals examples. These are unit tests, convention checks, and pedagogical demonstrations of operator algebra, perturbation theory, tensor conventions, symmetry, quadrature, and numerical differentiation.

## Numerical / algorithmic content

## Implementation structure

- Expansion relations for operator basis transforms.
- Accuracy threshold
- Test IST expansion of Zeeman level projectors
- Build known Zeeman projector
- Obtain IST expansion.
- Reconstruct operator from IST terms
- Report IST expansion failures
- Test BM expansion of bosonic level projectors
- Build known bosonic projector
- Obtain BM expansion.
- Reconstruct operator from BM terms
- Report BM expansion failures

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `irr_sph_ten()`, `enlev2ist()`, `coeffs()`, `states()`, `boson_mono()`, `enlev2bm()`, `weyl()`, `speye()`, `bos2ist()`, `sin_tran()`, `hdot()`, `oper2bm()`.
