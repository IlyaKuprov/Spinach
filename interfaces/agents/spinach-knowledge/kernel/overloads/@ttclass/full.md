# kernel/overloads/@ttclass/full.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/overloads/@ttclass/full.m`
- Signature: `answer=full(ttrain)`
- Total lines: 58

## Purpose

Converts a tensor train representation of a matrix into a matrix. Syntax: answer=full(ttrain)

## Physical / mathematical content

- Tensor-train linear algebra. These files implement compressed high-dimensional operators and AMEn/SVD-based algebra in tensor-train format.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 23-24: Preallocate the result; implemented by `answer=zeros(size(ttrain))`.
- Lines 26-27: Get object dimensions; implemented by `[ncores,ntrains]=size(ttrain.cores)`.
- Lines 29-30: Get tensor ranks; implemented by `ttranks=ranks(ttrain)`.
- Lines 32-33: Get mode sizes; implemented by `modesizes=sizes(ttrain)`.
- Lines 35-36: Loop over the buffer; implemented by `for n=1:ntrains`.
- Lines 38-39: Multiply up the tensor train; implemented by `next_term=reshape(ttrain.cores{ncores,n},[ttranks(ncores,n),modesizes(ncores,1)*modesizes(ncores,2)])`.
- Lines 47-48: Add the matrix to the result; implemented by `answer=answer+ttrain.coeff(n)*next_term`.

### Control flow inferred from the code

- Line 36: `for` loop over `n=1:ntrains`.
- Line 40: `for` loop over `k=(ncores-1):(-1):1`.

### Key state/data transformations

- Lines 24: computes `answer` using `answer=zeros(size(ttrain))`.
- Lines 27: computes `[ncores,ntrains]` using `[ncores,ntrains]=size(ttrain.cores)`.
- Lines 30: computes `ttranks` using `ttranks=ranks(ttrain)`.
- Lines 33: computes `modesizes` using `modesizes=sizes(ttrain)`.
- Lines 39: computes `next_term` using `next_term=reshape(ttrain.cores{ncores,n},[ttranks(ncores,n),modesizes(ncores,1)*modesizes(ncores,2)])`.

## Parameters / inputs

- ttrain -tensor train object

## Outputs

- answer -a full matrix
- Note: the result can be huge, careless use would crash the system.

## Implementation structure

- Converts a tensor train representation of a matrix
- into a matrix. Syntax:
- answer=full(ttrain)
- ttrain -tensor train object
- answer -a full matrix
- Note: the result can be huge, careless use would crash the system.
- Preallocate the result
- Get object dimensions
- Get tensor ranks
- Get mode sizes
- Loop over the buffer
- Multiply up the tensor train

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `ranks()`, `sizes()`, `ttranks()`, `modesizes()`.
