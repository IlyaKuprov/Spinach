# kernel/overloads/@ttclass/pack.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/overloads/@ttclass/pack.m`
- Signature: `ttout=pack(tt)`
- Total lines: 73

## Purpose

This subroutine packs all trains from the addition buffer into a single tensor train, but does not perform the recom- pression. Normally you should not call it directly, use ttclass/shrink.m instead. Syntax: ttout=pack(tt)

## Physical / mathematical content

- Tensor-train linear algebra. These files implement compressed high-dimensional operators and AMEn/SVD-based algebra in tensor-train format.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 24-25: Read tensor ranks and dimensions; implemented by `sz=tt.sizes`.
- Lines 29-30: Fast return if possible; implemented by `if N==1`.
- Lines 35-36: Total rank of all summands; implemented by `R=sum(rnk,2)`.
- Lines 39-40: Preallocate a single tensor train for all summands; implemented by `core=cell(d,1)`.
- Lines 45-46: Collect buffered trains in a single tensor train; implemented by `rr=[zeros(d+1,1), cumsum(rnk,2)]`.
- Lines 55-56: Special case of one-site tensor trains; implemented by `core{1,1}(1,:,:,1)=core{1,1}(1,:,:,1)+tt.coeff(1,n)*tt.cores{1,n}`.
- Lines 60-61: Construct empty ttclass array and fill the fields; implemented by `ttout=ttclass`.

### Control flow inferred from the code

- Line 30: conditional branch on `N==1`.
- Line 41: `for` loop over `k=1:d`.
- Line 47: `for` loop over `n=1:N`.
- Line 48: conditional branch on `d>1`.
- Line 50: `for` loop over `k=2:d-1`.

### Key state/data transformations

- Lines 25: computes `sz` using `sz=tt.sizes`.
- Lines 26: computes `rnk` using `rnk=tt.ranks`.
- Lines 27: computes `[d,N]` using `[d,N]=size(tt.cores)`.
- Lines 31: computes `ttout` using `ttout=tt`.
- Lines 36: computes `R` using `R=sum(rnk,2)`.
- Lines 37: computes `R(1)` using `R(1)=1; R(d+1)=1`.
- Lines 40: computes `core` using `core=cell(d,1)`.
- Lines 42: computes `core{k}` using `core{k}=zeros(R(k),sz(k,1),sz(k,2),R(k+1))`.
- Lines 46: computes `rr` using `rr=[zeros(d+1,1), cumsum(rnk,2)]`.
- Lines 49: computes `core{1,1}(1,:,:,rr(2,n)+1:rr(2,n+1))` using `core{1,1}(1,:,:,rr(2,n)+1:rr(2,n+1))=tt.coeff(1,n)*tt.cores{1,n}`.
- Lines 51: computes `core{k,1}(rr(k,n)+1:rr(k,n+1),:,:,rr(k+1,n)+1:rr(k+1,n+1))` using `core{k,1}(rr(k,n)+1:rr(k,n+1),:,:,rr(k+1,n)+1:rr(k+1,n+1))=tt.cores{k,n}`.
- Lines 53: computes `core{d,1}(rr(d,n)+1:rr(d,n+1),:,:,1)` using `core{d,1}(rr(d,n)+1:rr(d,n+1),:,:,1)=tt.cores{d,n}`.
- Lines 56: computes `core{1,1}(1,:,:,1)` using `core{1,1}(1,:,:,1)=core{1,1}(1,:,:,1)+tt.coeff(1,n)*tt.cores{1,n}`.
- Lines 62: computes `ttout.coeff` using `ttout.coeff=1`.
- Lines 63: computes `ttout.cores` using `ttout.cores=core`.
- Lines 64: computes `ttout.tolerance` using `ttout.tolerance=abs(sum(tt.tolerance))`.

## Parameters / inputs

- tt -tensor train object with unprocessed additions

## Outputs

- ttout -tensor train with additions buffer absorbed
- into the cores of the tensor, but not re-
- compressed

## Implementation structure

- This subroutine packs all trains from the addition buffer
- into a single tensor train, but does not perform the recom-
- pression. Normally you should not call it directly, use
- ttclass/shrink.m instead. Syntax:
- ttout=pack(tt)
- tt - tensor train object with unprocessed additions
- ttout -tensor train with additions buffer absorbed
- into the cores of the tensor, but not re-
- compressed
- Read tensor ranks and dimensions
- Fast return if possible
- Total rank of all summands

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `cumsum()`.
