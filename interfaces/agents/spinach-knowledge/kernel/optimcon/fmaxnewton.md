# kernel/optimcon/fmaxnewton.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/optimcon/fmaxnewton.m`
- Signature: `[x,data]=fmaxnewton(spin_system,cost_function,guess)`
- Total lines: 392

## Purpose

Finds a local maximum of a function of several variables using Newton and quasi-Newton algorithms. Syntax: [x,data]=fmaxnewton(spin_system,cost_function,guess)

## Physical / mathematical content

- Optimal-control core routines. These files implement GRAPE-style objective evaluation, quasi-Newton search, line search, regularisation, distortion models, and waveform parameterisations.
- The numerical method is quasi-Newton optimisation: curvature information is approximated from successive step and gradient differences instead of forming exact second derivatives every iteration.
- The numerical method is limited-memory quasi-Newton optimisation, which keeps only a short curvature history and is therefore suitable for waveform vectors too large for dense Hessians.
- The optimisation logic is Newton or Newton-like: search directions use first- and second-order local curvature information, usually with regularisation or line-search safeguards.

## Numerical / algorithmic content

- The code contains an inverse-problem or ill-conditioning aspect and therefore introduces explicit regularisation, model selection, or stabilisation logic.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `header()`, `footer()`, `itrep()`, `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 42-43: Check consistency; implemented by `grumble(spin_system,cost_function,guess)`.
- Lines 45-46: Initialise counters; implemented by `data.count.iter=0; data.count.fx=0`.
- Lines 50-53: Look for checkpoint file; implemented by `if isfield(spin_system.control,'checkpoint')&& exist([spin_system.sys.scratch filesep spin_system.control.checkpoint],'file')`.
- Lines 55-56: Read the initial guess from checkpoint file in global scratch; implemented by `report(spin_system,'WARNING: optimisation restarted from the checkpoint file')`.
- Lines 66-67: Stretch the guess supplied by user; implemented by `data.x_shape=size(guess); x=guess(:)`.
- Lines 71-72: Prepare the freeze mask; implemented by `if ~isempty(spin_system.control.freeze)`.
- Lines 74-75: Validate freeze mask; implemented by `if ~all(size(spin_system.control.freeze)==data.x_shape)`.
- Lines 79-80: Stretch freeze mask to match x; implemented by `frozen=spin_system.control.freeze(:)`.
- Lines 84-85: Do not freeze anything; implemented by `frozen=false(size(x))`.
- Lines 89-90: Print the header; implemented by `header(spin_system)`.
- Lines 92-93: Get the video writer going (use VLC to play); implemented by `if isfield(spin_system.control,'video_file')`.
- Lines 95-97: This format is available on all platforms; implemented by `VW=VideoWriter(spin_system.control.video_file, 'Motion JPEG AVI'); open(VW)`.
- Lines 101-103: If zero iterations, still display graphics, return fidelity and the procedural data; implemented by `if spin_system.control.max_iter==0`.
- Lines 107-108: Start the iteration loop; implemented by `for n=1:spin_system.control.max_iter`.
- Lines 110-111: Default exit flag; implemented by `exitflag=0`.
- Lines 113-114: Update iteration counter; implemented by `data.count.iter=data.count.iter+1`.
- Lines 116-117: Get the search direction; implemented by `switch spin_system.control.method`.
- Lines 123-124: Get objective and gradient; implemented by `[data,fx,g]=objeval(x,cost_function,data,spin_system)`.

### Control flow inferred from the code

- Line 51: conditional branch on `isfield(spin_system.control,'checkpoint')&&`.
- Line 59: conditional branch on `numel(x)~=numel(guess)`.
- Line 72: conditional branch on `~isempty(spin_system.control.freeze)`.
- Line 75: conditional branch on `~all(size(spin_system.control.freeze)==data.x_shape)`.
- Line 93: conditional branch on `isfield(spin_system.control,'video_file')`.
- Line 103: conditional branch on `spin_system.control.max_iter==0`.
- Line 108: `for` loop over `n=1:spin_system.control.max_iter`.
- Line 117: dispatches on `spin_system.control.method`; cases `{'lbfgs','rbfgs'}`, `{'newton','goodwin'}`.
- Line 121: conditional branch on `n==1`.
- Line 127: conditional branch on `(abs(data.fx_sep_pen(1))<1e-6)||(norm(g(~frozen),2)<1e-6)`.
- Line 145: conditional branch on `size(dx_hist,2)>spin_system.control.n_grads`.
- Line 148: conditional branch on `size(dg_hist,2)>spin_system.control.n_grads`.
- Line 154: conditional branch on `strcmp(spin_system.control.method,'lbfgs')`.
- Line 196: conditional branch on `any(~isfinite(dir))`.

### Key state/data transformations

- Lines 46: computes `data.count.iter` using `data.count.iter=0; data.count.fx=0`.
- Lines 47: computes `data.count.gfx` using `data.count.gfx=0; data.count.hfx=0`.
- Lines 48: computes `data.count.rfo` using `data.count.rfo=0`.
- Lines 62: computes `data.x_shape` using `data.x_shape=size(guess)`.
- Lines 80: computes `frozen` using `frozen=spin_system.control.freeze(:)`.
- Lines 96-97: computes `VW` using `VW=VideoWriter(spin_system.control.video_file, 'Motion JPEG AVI'); open(VW)`.
- Lines 104: computes `[data,fx]` using `[data,fx]=objeval(x,cost_function,data,spin_system)`.
- Lines 111: computes `exitflag` using `exitflag=0`.
- Lines 124: computes `[data,fx,g]` using `[data,fx,g]=objeval(x,cost_function,data,spin_system)`.
- Lines 132: computes `old_x` using `old_x=x(~frozen); dx_hist=[]`.
- Lines 133: computes `old_g` using `old_g=g(~frozen); dg_hist=[]`.
- Lines 136: computes `dir` using `dir=g.*(~frozen); dir=0.01*dir/max(abs(dir))`.
- Lines 141: computes `dx_hist` using `dx_hist=[x(~frozen)-old_x dx_hist]; old_x=x(~frozen)`.
- Lines 142: computes `dg_hist` using `dg_hist=[g(~frozen)-old_g dg_hist]; old_g=g(~frozen)`.
- Lines 157: computes `dir(~frozen)` using `dir(~frozen)=lbfgs(dx_hist,dg_hist,g(~frozen))`.
- Lines 162: computes `H` using `H=bfgs(dx_hist,dg_hist,g(~frozen))`.
- Lines 165: computes `[H,data]` using `[H,data]=hessreg(spin_system,H,g(~frozen),data)`.
- Lines 178: computes `[data,fx,g,H]` using `[data,fx,g,H]=objeval(x,cost_function,data,spin_system)`.

### Local helper functions

- Line 286: `header()` — `function header(spin_system)`. Footer printing function
  - Representative operation: `report(spin_system,'==============================================================================================')`.
  - Representative operation: `report(spin_system,'Iter #f #g #H #R fidelity penalties total alpha |grad| SDGA ')`.
- Line 293: `footer()` — `function footer(spin_system,exitflag,data)`.
  - Representative operation: `report(spin_system,'----------------------------------------------------------------------------------------------')`.
  - Representative operation: `switch(spin_system.control.method)`.
- Line 317: `itrep()` — `function itrep(spin_system,fx,g,dir,alpha,data)`. Performance figures
  - Representative operation: `fid=data.fx_sep_pen(1)`.
  - Representative operation: `pens=sum(data.fx_sep_pen(2:end))`.
- Line 342: `grumble()` — `function grumble(spin_system,cost_function,guess)`.
  - Representative operation: `if ~isa(cost_function,'function_handle')`.
  - Representative operation: `error('cost_function must be a function handle.')`.

## Parameters / inputs

- spin_system -Spinach data object that has been
- through optimcon.m function
- cost_function -a function handle that takes the input
- the size of guess
- guess -the initial point of the optimisation

## Outputs

- x -the final point of the optimisation
- data.count.iter -iteration counter
- data.count.fx -function evaluation counter
- data.count.gfx -gradient evaluation counter
- data.count.hfx -Hessian evaluation counter
- data.count.rfo -RFO iteration counter
- data.x_shape -output of size(guess)
- data.* -further fields may be set by the
- objective functon

## Implementation structure

- Finds a local maximum of a function of several variables using Newton
- and quasi-Newton algorithms. Syntax:
- [x,data]=fmaxnewton(spin_system,cost_function,guess)
- spin_system -Spinach data object that has been
- through optimcon.m function
- cost_function -a function handle that takes the input
- the size of guess
- guess -the initial point of the optimisation
- x -the final point of the optimisation
- data.count.iter -iteration counter
- data.count.fx -function evaluation counter
- data.count.gfx -gradient evaluation counter

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `isfield()`, `exist()`, `report()`, `load()`, `guess()`, `all()`, `int2str()`, `false()`, `header()`, `VideoWriter()`, `open()`, `objeval()`, `dx_hist()`, `dg_hist()`, `strcmp()`.
