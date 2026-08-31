# kernel/optimcon/bracketing.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/optimcon/bracketing.m`
- Signature: `[a,b,alpha,fx,gfx,next_act,data]=bracketing(cost_function,alpha,dir,x_0,...`
- Total lines: 177

## Purpose

Expands a trial step into a bracket that contains an acceptable line search point, or accepts the step directly if the Wolfe tests are met before sectioning becomes necessary.

## Physical / mathematical content

- Optimal-control core routines. These files implement GRAPE-style objective evaluation, quasi-Newton search, line search, regularisation, distortion models, and waveform parameterisations.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 55-56: Check consistency; implemented by `grumble(cost_function,alpha,dir,x_0,fx_0,gfx_0)`.
- Lines 58-59: Apply coordinate freezing mask when requested; implemented by `if ~isempty(spin_system.control.freeze)`.
- Lines 64-65: Evaluate objective and gradient at the first trial point; implemented by `[data,fx_2,gfx_2]=objeval(x_0+alpha*dir,cost_function,data,spin_system)`.
- Lines 72-73: Initialise empty bracket records; implemented by `a.alpha=[]; a.fx=[]; a.gfx=[]`.
- Lines 76-77: Initialise bracketing history variables; implemented by `fx=fx_0; fx_1=fx_0`.
- Lines 81-82: Expand bracket until acceptance or interval capture; implemented by `while true`.
- Lines 84-86: Capture bracket when Armijo or monotonicity fails; implemented by `if (~alpha_conds(1,alpha_2,fx_0,fx_2,gfx_0,[],dir,spin_system))|| (~alpha_conds(0,[],fx_0,fx_2,[],[],[],spin_system))`.
- Lines 88-89: Store current interval endpoints; implemented by `a.alpha=alpha_1; a.fx=fx_1; a.gfx=gfx_1`.
- Lines 92-93: Hand over to sectioning stage; implemented by `next_act='sectioning'; return`.
- Lines 97-98: Accept step immediately if curvature condition passes; implemented by `if alpha_conds(2,[],[],[],gfx_0,gfx_2,dir,spin_system)`.
- Lines 100-101: Return accepted step information; implemented by `alpha=alpha_2; fx=fx_2; gfx=gfx_2`.
- Lines 103-104: No further line search steps are required; implemented by `next_act='none'; return`.
- Lines 108-109: Capture bracket when directional derivative changes sign; implemented by `if ~alpha_conds(3,[],[],[],[],gfx_2,dir,spin_system)`.
- Lines 111-112: Store current interval endpoints; implemented by `a.alpha=alpha_2; a.fx=fx_2; a.gfx=gfx_2`.
- Lines 120-121: Build interpolation window ahead of current trial point; implemented by `br_end_pt_A=2*alpha_2-alpha_1`.
- Lines 124-125: Stop if the expansion has escaped to non-finite values; implemented by `if ~all(isfinite([br_end_pt_A br_end_pt_B fx_2 gfx_2'*dir]))`.
- Lines 129-131: Maximise the cubic model inside interpolation bounds; implemented by `alpha_new=cubic_interp(br_end_pt_A,br_end_pt_B,alpha_1,alpha_2, fx_1,gfx_1'*dir,fx_2,gfx_2'*dir)`.
- Lines 133-134: Shift history to the new trial point; implemented by `alpha_1=alpha_2; alpha_2=alpha_new`.

### Control flow inferred from the code

- Line 59: conditional branch on `~isempty(spin_system.control.freeze)`.
- Line 68: conditional branch on `~isempty(spin_system.control.freeze)`.
- Line 82: `while` loop over `true`.
- Line 85: conditional branch on `(~alpha_conds(1,alpha_2,fx_0,fx_2,gfx_0,[],dir,spin_system))||`.
- Line 98: conditional branch on `alpha_conds(2,[],[],[],gfx_0,gfx_2,dir,spin_system)`.
- Line 109: conditional branch on `~alpha_conds(3,[],[],[],[],gfx_2,dir,spin_system)`.
- Line 125: conditional branch on `~all(isfinite([br_end_pt_A br_end_pt_B fx_2 gfx_2'*dir]))`.
- Line 141: conditional branch on `~isempty(spin_system.control.freeze)`.

### Key state/data transformations

- Lines 60: computes `dir` using `dir=dir.*(~spin_system.control.freeze(:))`.
- Lines 61: computes `gfx_0` using `gfx_0=gfx_0.*(~spin_system.control.freeze(:))`.
- Lines 65: computes `[data,fx_2,gfx_2]` using `[data,fx_2,gfx_2]=objeval(x_0+alpha*dir,cost_function,data,spin_system)`.
- Lines 69: computes `gfx_2` using `gfx_2=gfx_2.*(~spin_system.control.freeze(:))`.
- Lines 73: computes `a.alpha` using `a.alpha=[]; a.fx=[]; a.gfx=[]`.
- Lines 74: computes `b.alpha` using `b.alpha=[]; b.fx=[]; b.gfx=[]`.
- Lines 77: computes `fx` using `fx=fx_0; fx_1=fx_0`.
- Lines 78: computes `gfx` using `gfx=gfx_0; gfx_1=gfx_0`.
- Lines 79: computes `alpha_1` using `alpha_1=0; alpha_2=alpha`.
- Lines 93: computes `next_act` using `next_act='sectioning'; return`.
- Lines 101: computes `alpha` using `alpha=alpha_2; fx=fx_2; gfx=gfx_2`.
- Lines 121: computes `br_end_pt_A` using `br_end_pt_A=2*alpha_2-alpha_1`.
- Lines 122: computes `br_end_pt_B` using `br_end_pt_B=alpha_2+spin_system.control.ls_tau1*(alpha_2-alpha_1)`.
- Lines 130-131: computes `alpha_new` using `alpha_new=cubic_interp(br_end_pt_A,br_end_pt_B,alpha_1,alpha_2, fx_1,gfx_1'*dir,fx_2,gfx_2'*dir)`.
- Lines 135: computes `fx_1` using `fx_1=fx_2; gfx_1=gfx_2; x_1=x_0+alpha_2*dir`.

### Local helper functions

- Line 150: `grumble()` — `function grumble(cost_function,alpha,dir,x_0,fx_0,gfx_0)`.
  - Representative operation: `if ~isa(cost_function,'function_handle')`.
  - Representative operation: `error('cost_function must be a function handle.')`.

## Syntax

```matlab
[A,B,alpha,fx,gfx,next_act,data]=...
bracketing(cost_function,alpha,dir,x_0,fx_0,...
gfx_0,data,spin_system)
```

## Parameters / inputs

- cost_function -objective function handle
- alpha -initial trial step length
- dir -search direction vector
- x_0 -current optimisation vector
- fx_0 -objective value at x_0
- gfx_0 -gradient at x_0
- data -optimisation workspace structure
- spin_system -Spinach data structure with
- line search settings

## Outputs

- A -lower bracket point structure
- B -upper bracket point structure
- alpha -accepted step length when found
- fx -objective value at accepted step
- gfx -gradient at accepted step
- next_act -continuation tag, either
- 'sectioning' or 'none'
- data -updated optimisation workspace

## Implementation structure

- Expands a trial step into a bracket that contains an acceptable
- line search point, or accepts the step directly if the Wolfe
- tests are met before sectioning becomes necessary.
- [A,B,alpha,fx,gfx,next_act,data]=...
- bracketing(cost_function,alpha,dir,x_0,fx_0,...
- gfx_0,data,spin_system)
- cost_function -objective function handle
- alpha -initial trial step length
- dir -search direction vector
- x_0 -current optimisation vector
- fx_0 -objective value at x_0
- gfx_0 -gradient at x_0

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `objeval()`, `alpha_conds()`, `all()`, `cubic_interp()`, `isscalar()`, `iscolumn()`, `isequal()`.
