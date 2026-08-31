# kernel/overloads/@ttclass/ttort.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/overloads/@ttclass/ttort.m`
- Signature: `[tt,lognrm]=ttort(tt,direct)`
- Total lines: 169

## Purpose

Performs TT-orthogonalisation for a tensor train (or for each tensor train in a buffered sum). Syntax: [tt,lognrm]=ttort(tt,direct)

## Physical / mathematical content

- Tensor-train linear algebra. These files implement compressed high-dimensional operators and AMEn/SVD-based algebra in tensor-train format.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 34-35: Read tensor ranks and dimensions; implemented by `sz=tt.sizes; rnk=tt.ranks`.
- Lines 38-39: Compute logs of norms; implemented by `if (nargout>1)`.
- Lines 44-45: Decide the direction; implemented by `switch direct`.
- Lines 47-48: Left-to-right direction; implemented by `case +1`.
- Lines 50-51: Loop over buffered trains; implemented by `for n=1:N`.
- Lines 53-54: Get the ranks of the current train; implemented by `r=rnk(:,n)`.
- Lines 56-57: TT-othogonalise left-to-right; implemented by `for k=1:d-1`.
- Lines 59-60: Shape cores into matrices; implemented by `B=reshape(tt.cores{k,n}, [r(k)*sz(k,1)*sz(k,2), r(k+1)])`.
- Lines 63-64: Run orthogonal-triangular decomposition; implemented by `[Q,R]=qr(B,0); nrm=norm(R,2)`.
- Lines 66-67: Take care of the norm; implemented by `if (nrm>0)`.
- Lines 76-77: Multiply triangular part forward; implemented by `C=R*C; rnew=size(Q,2)`.
- Lines 79-80: Fold the cores back; implemented by `tt.cores{k,n}=reshape(Q, [r(k), sz(k,1), sz(k,2), rnew])`.
- Lines 83-84: Update the ranks; implemented by `r(k+1)=rnew`.
- Lines 88-89: Take care of the last norm; implemented by `nrm=norm(tt.cores{d,n}(:),2)`.
- Lines 101-102: Right-to-left direction; implemented by `case -1`.
- Lines 110-111: TT-othogonalise right-to-left; implemented by `for k=d:-1:2`.
- Lines 113-114: Shape cores into matrices; implemented by `B=reshape(tt.cores{k,n}, [r(k), sz(k,1)*sz(k,2)*r(k+1)])`.
- Lines 117-118: Run orthogonal-triangular decomposition; implemented by `[Q,R]=qr(B.',0); nrm=norm(R,2)`.

### Control flow inferred from the code

- Line 39: conditional branch on `(nargout>1)`.
- Line 45: dispatches on `direct`; cases `+1`, `-1`.
- Line 51: `for` loop over `n=1:N`.
- Line 57: `for` loop over `k=1:d-1`.
- Line 67: conditional branch on `(nrm>0)`.
- Line 69: conditional branch on `(nargout>1)`.
- Line 90: conditional branch on `(nrm>0)`.
- Line 92: conditional branch on `(nargout>1)`.
- Line 105: `for` loop over `n=1:N`.
- Line 111: `for` loop over `k=d:-1:2`.
- Line 121: conditional branch on `(nrm>0)`.
- Line 123: conditional branch on `(nargout>1)`.
- Line 144: conditional branch on `(nrm>0)`.
- Line 146: conditional branch on `(nargout>1)`.

### Key state/data transformations

- Lines 35: computes `sz` using `sz=tt.sizes; rnk=tt.ranks`.
- Lines 36: computes `d` using `d=tt.ncores; N=tt.ntrains`.
- Lines 40: computes `lognrm` using `lognrm=log(tt.coeff)`.
- Lines 41: computes `tt.coeff` using `tt.coeff=ones(1,N)`.
- Lines 54: computes `r` using `r=rnk(:,n)`.
- Lines 60: computes `B` using `B=reshape(tt.cores{k,n}, [r(k)*sz(k,1)*sz(k,2), r(k+1)])`.
- Lines 61: computes `C` using `C=reshape(tt.cores{k+1,n}, [r(k+1), sz(k+1,1)*sz(k+1,2)*r(k+2)])`.
- Lines 64: computes `[Q,R]` using `[Q,R]=qr(B,0); nrm=norm(R,2)`.
- Lines 68: computes `R` using `R=R/nrm`.
- Lines 70: computes `lognrm(1,n)` using `lognrm(1,n)=lognrm(1,n)+log(nrm)`.
- Lines 72: computes `tt.coeff(1,n)` using `tt.coeff(1,n)=tt.coeff(1,n)*nrm`.
- Lines 80: computes `tt.cores{k,n}` using `tt.cores{k,n}=reshape(Q, [r(k), sz(k,1), sz(k,2), rnew])`.
- Lines 81: computes `tt.cores{k+1,n}` using `tt.cores{k+1,n}=reshape(C, [rnew, sz(k+1,1), sz(k+1,2), r(k+2)])`.
- Lines 84: computes `r(k+1)` using `r(k+1)=rnew`.
- Lines 89: computes `nrm` using `nrm=norm(tt.cores{d,n}(:),2)`.
- Lines 91: computes `tt.cores{d,n}` using `tt.cores{d,n}=tt.cores{d,n}/nrm`.
- Lines 135: computes `tt.cores{k-1,n}` using `tt.cores{k-1,n}=reshape(C, [r(k-1), sz(k-1,1), sz(k-1,2), rnew])`.
- Lines 138: computes `r(k)` using `r(k)=rnew`.

## Parameters / inputs

- direct=+1 -gives left-to-right orthogonality,
- direct=-1 -gives right-to-left orthogonality
- tt -tensor train object, possibly with buffered sums

## Outputs

- tt -tensor train object with all terms in the buffe-
- red sum has all of them orthogonalised in the
- direction requested
- lognrm -if this output is present, all buffered trains
- are also normalized, and natural logs of their
- norms returned in the vector lognrm. Use this
- option if the tensor norm is likely to exceed
- realmax()=1.7977e+308.
- Note: normally you should not call this subroutine directly.

## Implementation structure

- Performs TT-orthogonalisation for a tensor train (or for each tensor
- train in a buffered sum). Syntax:
- [tt,lognrm]=ttort(tt,direct)
- direct=+1 -gives left-to-right orthogonality,
- direct=-1 -gives right-to-left orthogonality
- tt -tensor train object, possibly with buffered sums
- tt -tensor train object with all terms in the buffe-
- red sum has all of them orthogonalised in the
- direction requested
- lognrm -if this output is present, all buffered trains
- are also normalized, and natural logs of their
- norms returned in the vector lognrm. Use this

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `rnk()`, `lognrm()`.
