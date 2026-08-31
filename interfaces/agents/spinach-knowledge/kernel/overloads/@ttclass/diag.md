# kernel/overloads/@ttclass/diag.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/overloads/@ttclass/diag.m`
- Signature: `tt=diag(tt)`
- Total lines: 73

## Purpose

Mimics the diag behaviour for tensor train matrix. Syntax: tt=diag(tt)

## Physical / mathematical content

- Tensor-train linear algebra. These files implement compressed high-dimensional operators and AMEn/SVD-based algebra in tensor-train format.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 22-23: Read tensor train sizes and ranks; implemented by `[ncores,ntrains]=size(tt.cores)`.
- Lines 26-27: Decide the dimensions; implemented by `if all(sz(:,1)==ones(ncores,1)) || all(sz(:,2)==ones(ncores,1))`.
- Lines 29-30: Vector on input, diagonal matrix on output; implemented by `for n=1:ntrains`.
- Lines 46-47: Matrix on input, column vector on output; implemented by `if all(sz(:,1)==sz(:,2))`.

### Control flow inferred from the code

- Line 27: conditional branch on `all(sz(:,1)==ones(ncores,1)) || all(sz(:,2)==ones(ncores,1))`.
- Line 30: `for` loop over `n=1:ntrains`.
- Line 31: `for` loop over `k=1:ncores`.
- Line 35: `for` loop over `p=1:r(k,n)`.
- Line 36: `for` loop over `q=1:r(k+1,n)`.
- Line 47: conditional branch on `all(sz(:,1)==sz(:,2))`.
- Line 48: `for` loop over `n=1:ntrains`.
- Line 49: `for` loop over `k=1:ncores`.
- Line 53: `for` loop over `p=1:r(k,n)`.
- Line 54: `for` loop over `q=1:r(k+1,n)`.

### Key state/data transformations

- Lines 23: computes `[ncores,ntrains]` using `[ncores,ntrains]=size(tt.cores)`.
- Lines 24: computes `r` using `r=ranks(tt); sz=sizes(tt)`.
- Lines 32: computes `vector_core` using `vector_core=tt.cores{k,n}`.
- Lines 33: computes `matrix_size` using `matrix_size=max(sz(k,1),sz(k,2))`.
- Lines 34: computes `matrix_core` using `matrix_core=zeros(r(k,n),matrix_size,matrix_size,r(k+1,n))`.
- Lines 37: computes `matrix_core(p,:,:,q)` using `matrix_core(p,:,:,q)=diag(reshape(vector_core(p,:,:,q),[sz(k,1),sz(k,2)]))`.
- Lines 40: computes `tt.cores{k,n}` using `tt.cores{k,n}=matrix_core`.
- Lines 51: computes `vector_size` using `vector_size=min(sz(k,1),sz(k,2))`.
- Lines 55: computes `vector_core(p,:,1,q)` using `vector_core(p,:,1,q)=diag(reshape(matrix_core(p,:,:,q),[sz(k,1),sz(k,2)]))`.

## Parameters / inputs

- tt -a tensor train representation of a matrix

## Outputs

- tt -if the input is a square matrix, returns a
- vector by computing diag of every core; if
- the input is a vector (one mode size is
- ones), returns a diagonal matrix

## Implementation structure

- Mimics the diag behaviour for tensor train matrix. Syntax:
- tt=diag(tt)
- tt -a tensor train representation of a matrix
- tt -if the input is a square matrix, returns a
- vector by computing diag of every core; if
- the input is a vector (one mode size is
- ones), returns a diagonal matrix
- Read tensor train sizes and ranks
- Decide the dimensions
- Vector on input, diagonal matrix on output
- Matrix on input, column vector on output
- The only mistake [the famous criminal finacier] Bernie Madoff

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `ranks()`, `sizes()`, `all()`, `matrix_core()`, `vector_core()`.
