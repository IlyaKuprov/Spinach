# kernel/overloads/@ttclass/ranks.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/overloads/@ttclass/ranks.m`
- Signature: `ttranks=ranks(ttrain)`
- Total lines: 45

## Purpose

Returns the bond dimensions of a tensor train. Syntax: ttranks=ranks(ttrain)

## Physical / mathematical content

- Tensor-train linear algebra. These files implement compressed high-dimensional operators and AMEn/SVD-based algebra in tensor-train format.
- The optimisation logic is Newton or Newton-like: search directions use first- and second-order local curvature information, usually with regularisation or line-search safeguards.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 22-23: Get core array dimensions; implemented by `[ncores,ntrains]=size(ttrain.cores)`.
- Lines 25-26: Preallocate the answer; implemented by `ttranks=zeros(ncores+1,ntrains)`.
- Lines 28-29: Loop over the buffer; implemented by `for n=1:ntrains`.
- Lines 31-32: Extract the ranks; implemented by `for k=1:ncores`.

### Control flow inferred from the code

- Line 29: `for` loop over `n=1:ntrains`.
- Line 32: `for` loop over `k=1:ncores`.

### Key state/data transformations

- Lines 23: computes `[ncores,ntrains]` using `[ncores,ntrains]=size(ttrain.cores)`.
- Lines 26: computes `ttranks` using `ttranks=zeros(ncores+1,ntrains)`.
- Lines 33: computes `ttranks(k,n)` using `ttranks(k,n)=size(ttrain.cores{k,n},1)`.
- Lines 35: computes `ttranks(ncores+1,n)` using `ttranks(ncores+1,n)=size(ttrain.cores{ncores,n},4)`.

## Parameters / inputs

- ttrain -a tensor train object

## Outputs

- ttranks -(ncores+1) by (ntrains) array; the
- first and the last elements for each
- train are 1

## Implementation structure

- Returns the bond dimensions of a tensor train. Syntax:
- ttranks=ranks(ttrain)
- ttrain -a tensor train object
- ttranks -(ncores+1) by (ntrains) array; the
- first and the last elements for each
- train are 1
- Get core array dimensions
- Preallocate the answer
- Loop over the buffer
- Extract the ranks
- I refrain from publishing for fear that disputes and controversies
- may be raised against me by ignoramuses.

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `ttranks()`.
