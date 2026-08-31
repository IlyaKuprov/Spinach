# kernel/overloads/@ttclass/amensum.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/overloads/@ttclass/amensum.m`
- Signature: `y=amensum(x,tol,opts)`
- Total lines: 388

## Purpose

Sums buffered tensor trains in a single tensor train using AMEn algorithm. Syntax: y=amensum(x,tol,opts)

## Physical / mathematical content

- Tensor-train linear algebra. These files implement compressed high-dimensional operators and AMEn/SVD-based algebra in tensor-train format.

## Numerical / algorithmic content

- The file also defines local helper function(s): `local_cp_to_tt()`, `step_tt_by_cp()`, `step_tt_by_tt()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 34-35: Process options; implemented by `if ~exist('opts','var'), opts = struct; end`.
- Lines 41-42: Read tensor train sizes and ranks; implemented by `[d,N]=size(x.cores); sz=x.sizes; rnk=x.ranks`.
- Lines 44-45: Proceed a CP format and a sum of TT formats separately; implemented by `if all(all(rnk==1))`.
- Lines 47-48: Convert ttclass to CP format; implemented by `X=cell(d,1)`.
- Lines 56-57: Generate a random initial guess; implemented by `[y,~]=ttort(rand(x, opts.init_guess_rank),-1)`.
- Lines 59-60: Generate a random enrichment if requested; implemented by `if opts.enrichment_rank>0`.
- Lines 64-65: Precompute interfaces of projections Y'X Z'X and Z'Y; implemented by `yx = cell(d+1,1); yx{1}=ones(1,N); yx{d+1}=ones(N,1)`.
- Lines 86-87: Set initial conditions; implemented by `flipped=false; satisfied=false; iter=0`.
- Lines 89-90: Main iteration cycle; implemented by `while ~satisfied`.
- Lines 92-93: Increment the counter; implemented by `iter=iter+1`.
- Lines 95-96: Read current ranks of approximation; implemented by `ry=y.ranks; sz=x.sizes`.
- Lines 101-102: The difference between two consequitive iterations; implemented by `largest_stepsize=0`.
- Lines 104-105: Loop over the cores; implemented by `for k=1:d`.
- Lines 107-108: Solve local approximation problem; implemented by `ynew=local_cp_to_tt(yx{k},yx{k+1},X{k,1},x.coeff)`.
- Lines 116-117: Check the convergence; implemented by `relative_stepsize=norm(ynew(:)-y.cores{k,1}(:),2)/norm(y.cores{k,1}(:),2)`.
- Lines 120-121: Run the truncation; implemented by `if k<d`.
- Lines 123-124: Truncate the updated core; implemented by `Y = reshape(ynew,[ry(k)*sz(k,1)*sz(k,2),ry(k+1)])`.
- Lines 132-133: Prepare error update and enrichment; implemented by `if opts.enrichment_rank>0`.

### Control flow inferred from the code

- Line 35: conditional branch on `~exist('opts','var'), opts = struct; end`.
- Line 36: conditional branch on `~isfield(opts, 'max_swp'); opts.max_swp=100; end`.
- Line 37: conditional branch on `~isfield(opts, 'init_guess_rank'); opts.init_guess_rank=2; end`.
- Line 38: conditional branch on `~isfield(opts, 'enrichment_rank'); opts.enrichment_rank=4; end`.
- Line 39: conditional branch on `~isfield(opts, 'verb'); opts.verb=0; end`.
- Line 45: conditional branch on `all(all(rnk==1))`.
- Line 49: `for` loop over `k=1:d`.
- Line 51: `for` loop over `i=1:N`.
- Line 60: conditional branch on `opts.enrichment_rank>0`.
- Line 66: conditional branch on `opts.enrichment_rank>0`.
- Line 71: `for` loop over `k=d:-1:2`.
- Line 73: conditional branch on `opts.enrichment_rank>0`.
- Line 80: conditional branch on `opts.enrichment_rank>0`.
- Line 90: `while` loop over `~satisfied`.

### Key state/data transformations

- Lines 42: computes `[d,N]` using `[d,N]=size(x.cores); sz=x.sizes; rnk=x.ranks`.
- Lines 48: computes `X` using `X=cell(d,1)`.
- Lines 50: computes `X{k,1}` using `X{k,1}=zeros(sz(k,1), sz(k,2),N)`.
- Lines 52: computes `X{k,1}(:,:,i)` using `X{k,1}(:,:,i)=reshape(x.cores{k,i}, [sz(k,1),sz(k,2),1])`.
- Lines 57: computes `[y,~]` using `[y,~]=ttort(rand(x, opts.init_guess_rank),-1)`.
- Lines 61: computes `[z,~]` using `[z,~]=ttort(rand(x, opts.enrichment_rank),-1)`.
- Lines 65: computes `yx` using `yx = cell(d+1,1); yx{1}=ones(1,N); yx{d+1}=ones(N,1)`.
- Lines 67: computes `zx` using `zx = cell(d+1,1); zx{1}=ones(1,N); zx{d+1}=ones(N,1)`.
- Lines 68: computes `zy` using `zy = cell(d+1,1); zy{1}=1; zy{d+1}=1`.
- Lines 70: computes `nrm` using `nrm = ones(d,1)`.
- Lines 72: computes `yx{k}` using `yx{k}=step_tt_by_cp(yx{k+1},y.cores{k,1},X{k,1},-1)`.
- Lines 74: computes `zx{k}` using `zx{k}=step_tt_by_cp(zx{k+1},z.cores{k,1},X{k,1},-1)`.
- Lines 75: computes `zy{k}` using `zy{k}=step_tt_by_tt(zy{k+1},z.cores{k,1},y.cores{k,1},-1)`.
- Lines 77: computes `nrm(k)` using `nrm(k)=norm(yx{k},'fro')`.
- Lines 87: computes `flipped` using `flipped=false; satisfied=false; iter=0`.
- Lines 93: computes `iter` using `iter=iter+1`.
- Lines 96: computes `ry` using `ry=y.ranks; sz=x.sizes`.
- Lines 98: computes `rz` using `rz=z.ranks`.

### Local helper functions

- Line 253: `local_cp_to_tt()` — `function ttcore = local_cp_to_tt(left_interface,right_interface,cp,coeff)`. Read sizes
  - Representative operation: `[r1,N1]=size(left_interface)`.
  - Representative operation: `[N2,r2]=size(right_interface)`.
- Line 285: `step_tt_by_cp()` — `function next_interface = step_tt_by_cp(interface, tt, cp, direction)`. Read sizes
  - Representative operation: `[r1,sz1,sz2,r2] = size(tt)`.
  - Representative operation: `[nn1,nn2,N] = size(cp)`.
- Line 337: `step_tt_by_tt()` — `function next_interface = step_tt_by_tt(interface, ttx, tty, direction)`. Read sizes
  - Representative operation: `[rx1,szx1,szx2,rx2] = size(ttx)`.
  - Representative operation: `[ry1,szy1,szy2,ry2] = size(tty)`.

## Parameters / inputs

- x -ttclass with buffered rank-one tensors
- tol -relative tolerance parameter, e.g. 1e-10
- The opts field is optional:
- opts.max_swp -maximum number of iterations
- opts.init_guess_rank -rank of the initial guess
- opts.enrichment_rank -rank of the enrichment
- opts.verb -verbosity switch

## Outputs

- y -ttclass with a single tensor train, such
- that |x-y|<tol*|x| in Frobenius norm

## Implementation structure

- Sums buffered tensor trains in a single tensor train using
- AMEn algorithm. Syntax:
- y=amensum(x,tol,opts)
- x -ttclass with buffered rank-one tensors
- tol -relative tolerance parameter, e.g. 1e-10
- The opts field is optional:
- opts.max_swp -maximum number of iterations
- opts.init_guess_rank -rank of the initial guess
- opts.enrichment_rank -rank of the enrichment
- opts.verb -verbosity switch
- y -ttclass with a single tensor train, such
- that |x-y|<tol*|x| in Frobenius norm

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `exist()`, `isfield()`, `all()`, `ttort()`, `step_tt_by_cp()`, `step_tt_by_tt()`, `nrm()`, `local_cp_to_tt()`, `ynew()`, `frob_chop()`, `znew()`, `enrichment()`, `ynew_enrich()`, `revert()`.
