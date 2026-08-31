# kernel/overloads/@ttclass/amensolve.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/overloads/@ttclass/amensolve.m`
- Signature: `x=amensolve(A,y,tol,opts,x0)`
- Total lines: 574

## Purpose

Solves the linear system Ax=y using the AMEn iteration. Syntax: x=amensolve(A,y,tol,opts,x0)

## Physical / mathematical content

- Tensor-train linear algebra. These files implement compressed high-dimensional operators and AMEn/SVD-based algebra in tensor-train format.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`, `revert_interface()`, `reduce_matrix()`, `reduce_vector()`, `local_matvec()`, `local_matrix()`, `local_vector()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 53-54: Get the default options started; implemented by `if (nargin<4)||isempty(opts), opts=struct; end`.
- Lines 56-57: Maximum number of AMEn sweeps; implemented by `if ~isfield(opts, 'nswp'); opts.nswp=20; end`.
- Lines 59-60: The rank of the initial guess; implemented by `if ~isfield(opts, 'init_guess_rank'); opts.init_guess_rank=2; end`.
- Lines 62-63: The rank of the residual and enrichment; implemented by `if ~isfield(opts, 'enrichment_rank'); opts.enrichment_rank=4; end`.
- Lines 65-66: Local accuracy gap; implemented by `if ~isfield(opts, 'resid_damp'); opts.resid_damp=2; end`.
- Lines 68-69: Maximum TT rank limit for the solution; implemented by `if ~isfield(opts, 'rmax'); opts.rmax=Inf; end`.
- Lines 71-72: Direct vs iterative solver switchover dimension; implemented by `if ~isfield(opts, 'max_full_size'); opts.max_full_size=500; end`.
- Lines 74-75: Maximum number of bicgstab iterations for local problems; implemented by `if ~isfield(opts, 'local_iters'); opts.local_iters=100; end`.
- Lines 77-78: Verbosity level: silent (0), sweep (1) or full (2); implemented by `if ~isfield(opts, 'verb'); opts.verb=1; end`.
- Lines 80-81: Initial guess; implemented by `if ~exist('x0','var')`.
- Lines 86-87: Check inputs for consistency; implemented by `grumble(A,y,x)`.
- Lines 89-90: Read dimension and mode sizes; implemented by `d=y.ncores`.
- Lines 94-95: Clear coefficients from all data; implemented by `A=clearcoeff(A)`.
- Lines 99-100: Initialize reductions of the system to the interfaces of the solution; implemented by `reduction_XAX = cell(d+1,1); reduction_XAX{1}=1; reduction_XAX{d+1}=1`.
- Lines 103-104: If the residual is used, initialise its ranks and reductions; implemented by `if opts.enrichment_rank>0`.
- Lines 106-107: Residual ranks; implemented by `rz = [1;opts.enrichment_rank*ones(d-1,1);1]`.
- Lines 109-110: Reductions of the system to the interfaces of the residual; implemented by `reduction_ZAX = cell(d+1,1); reduction_ZAX{1}=1; reduction_ZAX{d+1}=1`.
- Lines 115-117: The TT truncation accumulates the error from a single step with the factor sqrt(d). Therefore, the tolerance should be adjusted to sqrt(d).; implemented by `local_tolerance = tol/sqrt(d)`.

### Control flow inferred from the code

- Line 54: conditional branch on `(nargin<4)||isempty(opts), opts=struct; end`.
- Line 57: conditional branch on `~isfield(opts, 'nswp'); opts.nswp=20; end`.
- Line 60: conditional branch on `~isfield(opts, 'init_guess_rank'); opts.init_guess_rank=2; end`.
- Line 63: conditional branch on `~isfield(opts, 'enrichment_rank'); opts.enrichment_rank=4; end`.
- Line 66: conditional branch on `~isfield(opts, 'resid_damp'); opts.resid_damp=2; end`.
- Line 69: conditional branch on `~isfield(opts, 'rmax'); opts.rmax=Inf; end`.
- Line 72: conditional branch on `~isfield(opts, 'max_full_size'); opts.max_full_size=500; end`.
- Line 75: conditional branch on `~isfield(opts, 'local_iters'); opts.local_iters=100; end`.
- Line 78: conditional branch on `~isfield(opts, 'verb'); opts.verb=1; end`.
- Line 81: conditional branch on `~exist('x0','var')`.
- Line 104: conditional branch on `opts.enrichment_rank>0`.
- Line 129: `while` loop over `~satisfied`.
- Line 137: `for` loop over `k=d:-1:1`.
- Line 139: conditional branch on `~flipped`.

### Key state/data transformations

- Lines 82: computes `[x0,~]` using `[x0,~]=ttort(rand(y,opts.init_guess_rank),-1)`.
- Lines 84: computes `x` using `x=x0`.
- Lines 90: computes `d` using `d=y.ncores`.
- Lines 91: computes `sz` using `sz=A.sizes`.
- Lines 95: computes `A` using `A=clearcoeff(A)`.
- Lines 96: computes `y` using `y=clearcoeff(y)`.
- Lines 100: computes `reduction_XAX` using `reduction_XAX = cell(d+1,1); reduction_XAX{1}=1; reduction_XAX{d+1}=1`.
- Lines 101: computes `reduction_XY` using `reduction_XY = cell(d+1,1); reduction_XY{1}=1; reduction_XY{d+1}=1`.
- Lines 107: computes `rz` using `rz = [1;opts.enrichment_rank*ones(d-1,1);1]`.
- Lines 110: computes `reduction_ZAX` using `reduction_ZAX = cell(d+1,1); reduction_ZAX{1}=1; reduction_ZAX{d+1}=1`.
- Lines 111: computes `reduction_ZY` using `reduction_ZY = cell(d+1,1); reduction_ZY{1}=1; reduction_ZY{d+1}=1`.
- Lines 117: computes `local_tolerance` using `local_tolerance = tol/sqrt(d)`.
- Lines 120: computes `flipped` using `flipped = false`.
- Lines 121: computes `satisfied` using `satisfied = false`.
- Lines 122: computes `iter` using `iter = 1`.
- Lines 125: computes `norm_x` using `norm_x = ones(d-1,1)`.
- Lines 126: computes `norm_yAx` using `norm_yAx = ones(d-1,1)`.
- Lines 132: computes `ra` using `ra=A.ranks; ry=y.ranks; rx=x.ranks`.

### Local helper functions

- Line 412: `grumble()` — `function grumble(A,y,x0)`.
  - Representative operation: `if ~isa(A,'ttclass') || ~isa(y,'ttclass')`.
  - Representative operation: `error('Both matrix and right-hand-side should be ttclass.')`.
- Line 452: `revert_interface()` — `function [interface] = revert_interface(interface, ra, rw, rx)`.
  - Representative operation: `d = numel(interface)-1`.
  - Representative operation: `if (nargin<4)`.
- Line 469: `reduce_matrix()` — `function [interface,nrm] = reduce_matrix(interface, w, A, x)`.
  - Representative operation: `[ra1,n,m,ra2]=size(A)`.
  - Representative operation: `[rx1,~,~,rx2]=size(x)`.
- Line 498: `reduce_vector()` — `function [interface,nrm] = reduce_vector(interface, w, x)`.
  - Representative operation: `[rw1,n,~,rw2]=size(w)`.
  - Representative operation: `[rx1,~,~,rx2]=size(x)`.
- Line 514: `local_matvec()` — `function w=local_matvec(x, left_WAX, A, right_WAX)`.
  - Representative operation: `[ra1,n,m,ra2]=size(A)`.
  - Representative operation: `[rw1,rx1,~] = size(left_WAX)`.
- Line 535: `local_matrix()` — `function B=local_matrix(left_WAX, A, right_WAX)`.
  - Representative operation: `[ra1,n,m,ra2]=size(A)`.
  - Representative operation: `[rw1,rx1,~] = size(left_WAX)`.
- Line 553: `local_vector()` — `function w=local_vector(left_WX, x, right_WX)`. Had I the heavens' embroidered cloths,
  - Representative operation: `[rx1,n,~,rx2]=size(x)`.
  - Representative operation: `[rw1,~] = size(left_WX)`.

## Parameters / inputs

- A -ttclass representing a square matrix
- y -ttclass representing the right hand side
- tol -relative approximation and stopping tolerance,
- 1e-6 is a good start
- x0 -ttclass representing the initial guess
- Options (pass empty array for defaults):
- opts.nswp -maximum number of AMEn sweeps
- opts.init_guess_rank -the rank of the initial guess
- opts.enrichment_rank -the rank of the residual and
- enrichment
- opts.resid_damp -local accuracy gap
- opts.rmax -maximum TT rank limit for the solution
- opts.max_full_size -direct vs iterative solver switch-
- over dimension
- opts.local_iters -maximum number of bicgstab iterations
- for local problems
- opts.verb -Verbosity level: silent (0), sweep (1) or full (2)

## Outputs

- x -ttclass representing the solution such that
- |x-X| < tol |X| in Frobenius norm, where X is
- the exact solution.
- Note: A, y and x0 should have ntrains==1. Call shrink()
- on all three if that is not the case.

## Implementation structure

- Solves the linear system Ax=y using the AMEn iteration. Syntax:
- x=amensolve(A,y,tol,opts,x0)
- A -ttclass representing a square matrix
- y -ttclass representing the right hand side
- tol -relative approximation and stopping tolerance,
- 1e-6 is a good start
- x0 -ttclass representing the initial guess
- Options (pass empty array for defaults):
- opts.nswp -maximum number of AMEn sweeps
- opts.init_guess_rank -the rank of the initial guess
- opts.enrichment_rank -the rank of the residual and
- enrichment

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `isfield()`, `exist()`, `ttort()`, `grumble()`, `clearcoeff()`, `local_vector()`, `local_matvec()`, `local_matrix()`, `local_iters()`, `bicgstab()`, `frob_chop()`, `norm_x()`, `reduce_matrix()`, `reduce_vector()`, `norm_yAx()`, `revert_interface()`.
