# kernel/steady.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/steady.m`
- Signature: `rho=steady(spin_system,P,rho,method)`
- Total lines: 220

## Purpose

Steady state under the repeated action by the same dissi- pative evolution propagator. Syntax: rho=steady(spin_system,P,rho,tol,method)

## Physical / mathematical content

- This file belongs to the `kernel` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.
- The optimisation logic is Newton or Newton-like: search directions use first- and second-order local curvature information, usually with regularisation or line-search safeguards.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `size()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 40-41: Default initial guess; implemented by `if (~exist('rho','var'))||isempty(rho)`.
- Lines 51-52: Default method; implemented by `if (~exist('method','var'))||isempty(method)`.
- Lines 56-57: Check consistency; implemented by `grumble(spin_system,P,rho,method)`.
- Lines 59-60: Pick the method; implemented by `switch method`.
- Lines 64-65: Iteration stats; implemented by `cheap_diff=1; n_iter=1`.
- Lines 67-68: Trace functional of the stretched density matrix; implemented by `if strcmp(spin_system.bas.formalism,'zeeman-liouv')`.
- Lines 72-73: Keep going; implemented by `while cheap_diff>spin_system.tols.stst_tol`.
- Lines 75-76: Compute the square; implemented by `Psq=clean_up(spin_system,P*P,spin_system.tols.prop_chop)`.
- Lines 78-79: Pin the trace conservation row exactly; implemented by `if strcmp(spin_system.bas.formalism,'zeeman-liouv')`.
- Lines 83-84: Compute the difference and close the loop; implemented by `cheap_diff=max(abs(Psq-P),[],'all')`.
- Lines 87-88: Detect algorithm stagnation; implemented by `if n_iter>30, error('steady state convergence failure.'); end`.
- Lines 92-93: Compute the state; implemented by `rho=P*rho`.
- Lines 97-98: Get the Jacobian; implemented by `J=P-speye(size(P)); du=1`.
- Lines 100-101: Pick the normalisation pinning strategy; implemented by `switch spin_system.bas.formalism`.
- Lines 105-106: Pre-factor the Jacobian; implemented by `[LF,UF,RP]=lu(J(2:end,2:end))`.
- Lines 108-109: Iteration counter; implemented by `n_iter=0`.
- Lines 111-112: Newton iteration with unit state pinning; implemented by `while norm(du,2)>spin_system.tols.stst_tol`.
- Lines 114-115: Compute the residual; implemented by `r=P*rho-rho; r=r(2:end)`.

### Control flow inferred from the code

- Line 41: conditional branch on `(~exist('rho','var'))||isempty(rho)`.
- Line 42: dispatches on `spin_system.bas.formalism`; cases `'sphten-liouv'`, `'zeeman-liouv'`, `'squaring'`, `'newton'`, `'sphten-liouv'`, `'zeeman-liouv'`.
- Line 52: conditional branch on `(~exist('method','var'))||isempty(method)`.
- Line 60: dispatches on `method`; cases `'squaring'`, `'newton'`, `'sphten-liouv'`, `'zeeman-liouv'`.
- Line 68: conditional branch on `strcmp(spin_system.bas.formalism,'zeeman-liouv')`.
- Line 73: `while` loop over `cheap_diff>spin_system.tols.stst_tol`.
- Line 79: conditional branch on `strcmp(spin_system.bas.formalism,'zeeman-liouv')`.
- Line 88: conditional branch on `n_iter>30, error('steady state convergence failure.'); end`.
- Line 101: dispatches on `spin_system.bas.formalism`; cases `'sphten-liouv'`, `'zeeman-liouv'`.
- Line 112: `while` loop over `norm(du,2)>spin_system.tols.stst_tol`.
- Line 124: conditional branch on `n_iter>10, error('steady state convergence failure.'); end`.
- Line 143: `while` loop over `norm(du,2)>spin_system.tols.stst_tol`.
- Line 155: conditional branch on `n_iter>10, error('steady state convergence failure.'); end`.

### Key state/data transformations

- Lines 44: computes `rho` using `rho=zeros([size(P,2) 1],'like',1i); rho(1)=1`.
- Lines 46: computes `dim` using `dim=sqrt(size(P,2))`.
- Lines 53: computes `method` using `method='newton'`.
- Lines 65: computes `cheap_diff` using `cheap_diff=1; n_iter=1`.
- Lines 76: computes `Psq` using `Psq=clean_up(spin_system,P*P,spin_system.tols.prop_chop)`.
- Lines 85: computes `P` using `P=Psq; n_iter=n_iter+1`.
- Lines 98: computes `J` using `J=P-speye(size(P)); du=1`.
- Lines 106: computes `[LF,UF,RP]` using `[LF,UF,RP]=lu(J(2:end,2:end))`.
- Lines 109: computes `n_iter` using `n_iter=0`.
- Lines 115: computes `r` using `r=P*rho-rho; r=r(2:end)`.
- Lines 118: computes `du` using `du=-UF\(LF\(RP*r))`.
- Lines 121: computes `rho(2:end)` using `rho(2:end)=rho(2:end)+du; n_iter=n_iter+1`.
- Lines 134: computes `B` using `B=[J u0; u0' 0]`.

### Local helper functions

- Line 171: `grumble()` — `function grumble(spin_system,P,rho,method)`.
  - Representative operation: `if ~ismember(spin_system.bas.formalism,{'sphten-liouv','zeeman-liouv'})`.
  - Representative operation: `error('steady state is only available for sphten-liouv and zeeman-liouv formalisms.')`.

## Parameters / inputs

- P -propagator, an exponential of the Liouvillian
- that contains a thermalised relaxation super-
- operator (inter.equilibrium='IME' or 'dibari')
- or a product thereof (for example, from a re-
- peating block of a pulse or a pulse sequence)
- rho -optional initial guess for the steady state,
- a good one can significantly accelerate this
- function (leave empty otherwise); the state
- must have unit trace, which in sphten-liouv
- means a first element equal to 1
- method -'newton' (default) for the Newton-Raphson
- steady state solver, 'squaring' for propa-
- gator squaring (much more expensive, but
- unconditionally numerically stable)

## Outputs

- rho -steady state under the repeated applicati-
- on of the propagator P
- Note: available for sphten-liouv and zeeman-liouv formalisms; the
- Newton-Raphson solver pins the first state vector element in
- sphten-liouv and the density matrix trace in zeeman-liouv.

## Implementation structure

- Steady state under the repeated action by the same dissi-
- pative evolution propagator. Syntax:
- rho=steady(spin_system,P,rho,tol,method)
- P -propagator, an exponential of the Liouvillian
- that contains a thermalised relaxation super-
- operator (inter.equilibrium='IME' or 'dibari')
- or a product thereof (for example, from a re-
- peating block of a pulse or a pulse sequence)
- rho -optional initial guess for the steady state,
- a good one can significantly accelerate this
- function (leave empty otherwise); the state
- must have unit trace, which in sphten-liouv

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `exist()`, `rho()`, `speye()`, `complex()`, `grumble()`, `strcmp()`, `clean_up()`, `ismember()`, `ischar()`, `iscolumn()`.
