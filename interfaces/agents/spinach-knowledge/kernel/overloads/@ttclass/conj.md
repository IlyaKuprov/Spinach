# kernel/overloads/@ttclass/conj.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/overloads/@ttclass/conj.m`
- Signature: `tt=conj(tt)`
- Total lines: 44

## Purpose

Conjugates all core elements and coefficients of a tensor train object. Syntax: tt=conj(tt)

## Physical / mathematical content

- Tensor-train linear algebra. These files implement compressed high-dimensional operators and AMEn/SVD-based algebra in tensor-train format.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 22-23: Read tensor train sizes and ranks; implemented by `[ncores,ntrains]=size(tt.cores)`.
- Lines 25-26: Conjugate the cores; implemented by `for n=1:ntrains`.
- Lines 32-33: Conjugate the coefficients; implemented by `tt.coeff=conj(tt.coeff)`.

### Control flow inferred from the code

- Line 26: `for` loop over `n=1:ntrains`.
- Line 27: `for` loop over `k=1:ncores`.

### Key state/data transformations

- Lines 23: computes `[ncores,ntrains]` using `[ncores,ntrains]=size(tt.cores)`.
- Lines 28: computes `tt.cores{k,n}` using `tt.cores{k,n}=conj(tt.cores{k,n})`.
- Lines 33: computes `tt.coeff` using `tt.coeff=conj(tt.coeff)`.

## Parameters / inputs

- tt -tensor train object

## Outputs

- tt -tensor train object with complex-conjugated cores
- and coefficients

## Implementation structure

- Conjugates all core elements and coefficients of a tensor
- train object. Syntax:
- tt=conj(tt)
- tt -tensor train object
- tt -tensor train object with complex-conjugated cores
- and coefficients
- Read tensor train sizes and ranks
- Conjugate the cores
- Conjugate the coefficients
- I asked God for a bike, but I know God
- doesn't work that way. So I stole a bike
- and asked for forgiveness.
