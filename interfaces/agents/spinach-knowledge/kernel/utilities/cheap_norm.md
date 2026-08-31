# kernel/utilities/cheap_norm.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/cheap_norm.m`
- Signature: `n=cheap_norm(A,t,itmax)`
- Total lines: 182

## Purpose

The cheapest norm for various representations of matrices. CUDA stores matrices by rows, Matlab by columns, and polyadic objects can only multiply vectors, in which case Algorithm 2.4 from Hig- ham and Tisseur's paper: is used. Syntax: n=cheap_norm(A,t,itmax)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 36-37: Set the defaults; implemented by `if nargin<2, t=1; end`.
- Lines 40-41: Check consistency; implemented by `grumble(A,t,itmax)`.
- Lines 43-44: Inf-norm is best on GPU; implemented by `if isa(A,'gpuArray')`.
- Lines 48-49: 1-norm is best on CPU; implemented by `if ~isa(A,'polyadic')`.
- Lines 53-54: Relevant dimensions; implemented by `row_dim=size(A,1); col_dim=size(A,2)`.
- Lines 56-57: Respect small dimensions; implemented by `t=min(t,col_dim)`.
- Lines 59-60: Starting matrix; implemented by `X=ones(col_dim,t)`.
- Lines 80-81: Estimator state; implemented by `idx_hist=[]; idx=1:t; idx_best=1; est_old=0; S=zeros(row_dim,t)`.
- Lines 83-84: Iteration loop; implemented by `for k=1:(itmax+1)`.
- Lines 86-87: Matrix-vector products with the operator; implemented by `Y=A*X; col_norms=sum(abs(Y),1)`.
- Lines 89-90: Extract the current estimate; implemented by `[est,best_col]=max(col_norms)`.
- Lines 95-96: Stop if the estimate has not improved; implemented by `if (k>=2)&&(est<=est_old)`.
- Lines 100-101: Keep track of the best estimate; implemented by `est_old=est; S_old=S`.
- Lines 106-107: Phase matrix with the zero convention from the paper; implemented by `S=sign(Y); S(S==0)=1`.
- Lines 109-110: Remove redundant sign probes in the real case; implemented by `if isreal(A)`.
- Lines 130-131: Adjoint products with the phase matrix; implemented by `Z=A'*S; row_scores=max(abs(Z),[],2)`.
- Lines 133-134: Stop when the best column has been reached; implemented by `if (k>=2)&&(max(row_scores)==row_scores(idx_best))`.
- Lines 138-139: Follow the largest row scores; implemented by `[~,idx]=sort(row_scores,'descend'); idx=idx(:).'`.

### Control flow inferred from the code

- Line 37: conditional branch on `nargin<2, t=1; end`.
- Line 38: conditional branch on `nargin<3, itmax=5; end`.
- Line 44: conditional branch on `isa(A,'gpuArray')`.
- Line 49: conditional branch on `~isa(A,'polyadic')`.
- Line 61: conditional branch on `t>1`.
- Line 62: `for` loop over `col_idx=2:t`.
- Line 65: `for` loop over `col_idx=1:t`.
- Line 68: `while` loop over `any(abs(sign_pool'*X(:,col_idx))==col_dim)`.
- Line 72: conditional branch on `ntries>100`.
- Line 84: `for` loop over `k=1:(itmax+1)`.
- Line 91: conditional branch on `(est>est_old)||(k==2)`.
- Line 92: conditional branch on `k>=2, idx_best=idx(best_col); end`.
- Line 96: conditional branch on `(k>=2)&&(est<=est_old)`.
- Line 102: conditional branch on `k>itmax`.

### Key state/data transformations

- Lines 45: computes `n` using `n=norm(A,inf); return`.
- Lines 54: computes `row_dim` using `row_dim=size(A,1); col_dim=size(A,2)`.
- Lines 57: computes `t` using `t=min(t,col_dim)`.
- Lines 60: computes `X` using `X=ones(col_dim,t)`.
- Lines 63: computes `X(:,col_idx)` using `X(:,col_idx)=2*randi([0 1],col_dim,1)-1`.
- Lines 66: computes `sign_pool` using `sign_pool=X(:,[1:(col_idx-1) (col_idx+1):t])`.
- Lines 67: computes `ntries` using `ntries=0`.
- Lines 81: computes `idx_hist` using `idx_hist=[]; idx=1:t; idx_best=1; est_old=0; S=zeros(row_dim,t)`.
- Lines 87: computes `Y` using `Y=A*X; col_norms=sum(abs(Y),1)`.
- Lines 90: computes `[est,best_col]` using `[est,best_col]=max(col_norms)`.
- Lines 101: computes `est_old` using `est_old=est; S_old=S`.
- Lines 119: computes `S(:,col_idx)` using `S(:,col_idx)=2*randi([0 1],row_dim,1)-1`.
- Lines 131: computes `Z` using `Z=A'*S; row_scores=max(abs(Z),[],2)`.
- Lines 139: computes `[~,idx]` using `[~,idx]=sort(row_scores,'descend'); idx=idx(:).'`.
- Lines 146: computes `idx` using `idx=[idx(~ismember(idx,idx_hist)) idx(ismember(idx,idx_hist))]`.
- Lines 152: computes `X(idx(col_idx),col_idx)` using `X(idx(col_idx),col_idx)=1`.

### Local helper functions

- Line 164: `grumble()` — `function grumble(A,t,itmax)`.
  - Representative operation: `if ~isnumeric(A)`.
  - Representative operation: `error('A must be a matrix or a polyadic.')`.

## Parameters / inputs

- A -a matrix, or a polyadic representation thereof
- t -(optional) number of probe columns in the poly-
- adic norm estimator, defaults to 1
- itmax -(optional) maximum number of estimator iterati-
- ons, defaults to 5

## Outputs

- n -infinity-norm for GPU arrays, 1-norm for CPU arrays,
- and a lower-bound 1-norm estimate for polyadics
- Note: some norms are vastly more expensive than others, this
- function uses the cheapest ones available.

## Implementation structure

- The cheapest norm for various representations of matrices. CUDA
- stores matrices by rows, Matlab by columns, and polyadic objects
- can only multiply vectors, in which case Algorithm 2.4 from Hig-
- ham and Tisseur's paper:
- is used. Syntax:
- n=cheap_norm(A,t,itmax)
- A -a matrix, or a polyadic representation thereof
- t -(optional) number of probe columns in the poly-
- adic norm estimator, defaults to 1
- itmax -(optional) maximum number of estimator iterati-
- ons, defaults to 5
- n -infinity-norm for GPU arrays, 1-norm for CPU arrays,

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `randi()`, `any()`, `idx()`, `sign()`, `all()`, `row_scores()`, `ismember()`, `isscalar()`.
