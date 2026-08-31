# examples/fundamentals/convention_tests/tensors.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/fundamentals/convention_tests/tensors.m`
- Signature: `tensors()`
- Total lines: 45

## Purpose

Test the conversion from Stevens operator coefficients to irreducible spherical tensor coefficients.

## Physical / mathematical content

- Fundamentals examples. These are unit tests, convention checks, and pedagogical demonstrations of operator algebra, perturbation theory, tensor conventions, symmetry, quadrature, and numerical differentiation.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 8-9: Test up to rank 6 on spin 15/2; implemented by `max_rank=6; mult=16`.
- Lines 11-12: Generate a random set of Stevens operator coefficients; implemented by `for k=1:max_rank, r{k}=rand(2*k+1,1)-0.5; end`.
- Lines 14-15: Build the linear combination; implemented by `lin_comb_a=sparse(0)`.
- Lines 22-23: Translate the coefficients into ISTs; implemented by `for k=1:max_rank, r{k}=stev2sph(k,r{k}); end`.
- Lines 25-26: Build the linear combination; implemented by `lin_comb_b=sparse(0)`.
- Lines 34-35: Subtract the matrices and check the norm; implemented by `if norm(lin_comb_b-lin_comb_a,1)<1e-6`.

### Control flow inferred from the code

- Line 12: `for` loop over `k=1:max_rank, r{k}=rand(2*k+1,1)-0.5; end`.
- Line 16: `for` loop over `k=1:max_rank`.
- Line 17: `for` loop over `q=1:(2*k+1)`.
- Line 23: `for` loop over `k=1:max_rank, r{k}=stev2sph(k,r{k}); end`.
- Line 27: `for` loop over `k=1:max_rank`.
- Line 29: `for` loop over `q=1:(2*k+1)`.
- Line 35: conditional branch on `norm(lin_comb_b-lin_comb_a,1)<1e-6`.

### Key state/data transformations

- Lines 9: computes `max_rank` using `max_rank=6; mult=16`.
- Lines 15: computes `lin_comb_a` using `lin_comb_a=sparse(0)`.
- Lines 26: computes `lin_comb_b` using `lin_comb_b=sparse(0)`.
- Lines 28: computes `T` using `T=irr_sph_ten(mult,k)`.

## Implementation structure

- Test the conversion from Stevens operator coefficients to
- irreducible spherical tensor coefficients.
- Test up to rank 6 on spin 15/2
- Generate a random set of Stevens operator coefficients
- Build the linear combination
- Translate the coefficients into ISTs
- Subtract the matrices and check the norm

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `stevens()`, `stev2sph()`, `irr_sph_ten()`.
