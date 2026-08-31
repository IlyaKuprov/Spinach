# kernel/overloads/@ttclass/truncate.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/overloads/@ttclass/truncate.m`
- Signature: `ttout=truncate(tt)`
- Total lines: 72

## Purpose

Performs right-to-left SVD recompression for a tensor train. This should not be called directly, use shrink.m instead. Syntax: ttout=truncate(tt)

## Physical / mathematical content

- Tensor-train linear algebra. These files implement compressed high-dimensional operators and AMEn/SVD-based algebra in tensor-train format.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 25-26: Check consistency; implemented by `grumble(tt)`.
- Lines 28-29: Read tensor ranks and dimensions; implemented by `sz=tt.sizes; rnk=tt.ranks; d=tt.ncores`.
- Lines 31-32: Preallocate the result; implemented by `ttout=tt; r=rnk(:,1)`.
- Lines 34-35: Define the relative approximation accuracy; implemented by `eps=tt.tolerance(1,1)/tt.coeff(1,1)`.
- Lines 38-39: SVD cores right-to-left; implemented by `for k=d:-1:2`.
- Lines 50-51: Normalise the output; implemented by `nrm=norm(ttout.cores{1,1}(:),2)`.

### Control flow inferred from the code

- Line 39: `for` loop over `k=d:-1:2`.

### Key state/data transformations

- Lines 29: computes `sz` using `sz=tt.sizes; rnk=tt.ranks; d=tt.ncores`.
- Lines 32: computes `ttout` using `ttout=tt; r=rnk(:,1)`.
- Lines 35: computes `eps` using `eps=tt.tolerance(1,1)/tt.coeff(1,1)`.
- Lines 40: computes `C` using `C=reshape(ttout.cores{k,1}, [r(k), sz(k,1)*sz(k,2)*r(k+1)])`.
- Lines 41: computes `B` using `B=reshape(ttout.cores{k-1,1}, [r(k-1)*sz(k-1,1)*sz(k-1,2), r(k)])`.
- Lines 42: computes `[U,S,V]` using `[U,S,V]=svd(C,'econ'); S=diag(S)`.
- Lines 43: computes `rnew` using `rnew=frob_chop(S,eps*norm(S,'fro'))`.
- Lines 44: computes `U` using `U=U(:,1:rnew); S=S(1:rnew); V=V(:,1:rnew)`.
- Lines 45: computes `ttout.cores{k,1}` using `ttout.cores{k,1}=reshape(V', [rnew, sz(k,1), sz(k,2), r(k+1)])`.
- Lines 46: computes `ttout.cores{k-1,1}` using `ttout.cores{k-1,1}=reshape(B*U*diag(S), [r(k-1), sz(k-1,1), sz(k-1,2), rnew])`.
- Lines 47: computes `r(k)` using `r(k)=rnew`.
- Lines 51: computes `nrm` using `nrm=norm(ttout.cores{1,1}(:),2)`.
- Lines 52: computes `ttout.cores{1,1}` using `ttout.cores{1,1}=ttout.cores{1,1}/nrm`.
- Lines 53: computes `ttout.coeff(1,1)` using `ttout.coeff(1,1)=ttout.coeff(1,1)*nrm`.

### Local helper functions

- Line 58: `grumble()` — `function grumble(tt)`. "Even with basic 128-bit hash functions, such as MD5, [a hash collision] is a vanishingly rare event in the physical scien-
  - Representative operation: `if tt.ntrains>1`.
  - Representative operation: `error('Please sum all tensor trains before calling truncate')`.

## Parameters / inputs

- tt -a tensor train object with tt.ntrains=1
- and orthogonalised left-to-right

## Outputs

- ttout -a tensor train object, orthogonalised
- right-to-left
- Note: approximation tolerance (in Frobenius norm) is read from
- tt.tolerance property.

## Implementation structure

- Performs right-to-left SVD recompression for a tensor train. This
- should not be called directly, use shrink.m instead. Syntax:
- ttout=truncate(tt)
- tt -a tensor train object with tt.ntrains=1
- and orthogonalised left-to-right
- ttout -a tensor train object, orthogonalised
- right-to-left
- Note: approximation tolerance (in Frobenius norm) is read from
- tt.tolerance property.
- Check consistency
- Read tensor ranks and dimensions
- Preallocate the result

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `rnk()`, `frob_chop()`.
