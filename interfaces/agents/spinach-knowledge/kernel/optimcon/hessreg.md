# kernel/optimcon/hessreg.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/optimcon/hessreg.m`
- Signature: `[H,data]=hessreg(spin_system,H,g,data)`
- Total lines: 101

## Purpose

RFO regularisation for Newton-Raphson Hessian and gradient pairs. Syntax: [H,data]=hessreg(spin_system,H,g,data)

## Physical / mathematical content

- Optimal-control core routines. These files implement GRAPE-style objective evaluation, quasi-Newton search, line search, regularisation, distortion models, and waveform parameterisations.
- The optimisation logic is Newton or Newton-like: search directions use first- and second-order local curvature information, usually with regularisation or line-search safeguards.

## Numerical / algorithmic content

- An eigenvalue problem is solved or analysed, so the file is extracting spectra, stationary states, avoided crossings, or modal structure from the effective Hamiltonian or superoperator.
- The code contains an inverse-problem or ill-conditioning aspect and therefore introduces explicit regularisation, model selection, or stabilisation logic.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `size()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 29-30: Check consistency; implemented by `grumble(H,g)`.
- Lines 32-33: Set shorthands; implemented by `alpha=spin_system.control.reg_alpha`.
- Lines 38-39: Check if RFO is needed; implemented by `[~,p]=chol(H); pos_def=~logical(p)`.
- Lines 44-45: Start RFO iteration loop; implemented by `for ind=1:max_iter`.
- Lines 47-48: Make auxiliary Hessian; implemented by `H=[(alpha^2)*H alpha*g`.
- Lines 51-52: Calculate eigenvalue shift; implemented by `sigma=min([0 min(eig(H,'vector'))])`.
- Lines 54-55: Apply eigenvalue shift; implemented by `H=H-sigma*speye(size(H))`.
- Lines 57-58: Truncate and scale back; implemented by `H=H(1:(end-1),(1:end-1))`.
- Lines 61-62: Update alpha; implemented by `alpha=alpha*phi`.
- Lines 64-65: Increment the counter; implemented by `data.count.rfo=data.count.rfo+1`.
- Lines 67-68: Break if condition number reached; implemented by `if cond(H,2)<max_cond, break; end`.
- Lines 72-73: Clean up the result; implemented by `H=real(H+H')/2`.
- Lines 75-76: Warn if the target condition number was not reached; implemented by `if cond(H,2)>=max_cond`.

### Control flow inferred from the code

- Line 40: conditional branch on `pos_def&&(cond(H,2)<max_cond)`.
- Line 45: `for` loop over `ind=1:max_iter`.
- Line 68: conditional branch on `cond(H,2)<max_cond, break; end`.
- Line 76: conditional branch on `cond(H,2)>=max_cond`.

### Key state/data transformations

- Lines 33: computes `alpha` using `alpha=spin_system.control.reg_alpha`.
- Lines 34: computes `phi` using `phi=spin_system.control.reg_phi`.
- Lines 35: computes `max_iter` using `max_iter=spin_system.control.reg_max_iter`.
- Lines 36: computes `max_cond` using `max_cond=spin_system.control.reg_max_cond`.
- Lines 39: computes `[~,p]` using `[~,p]=chol(H); pos_def=~logical(p)`.
- Lines 48: computes `H` using `H=[(alpha^2)*H alpha*g`.
- Lines 52: computes `sigma` using `sigma=min([0 min(eig(H,'vector'))])`.
- Lines 65: computes `data.count.rfo` using `data.count.rfo=data.count.rfo+1`.

### Local helper functions

- Line 83: `grumble()` — `function grumble(H,g)`.
  - Representative operation: `if (~isnumeric(H))||(size(H,1)~=size(H,2))|| (~issymmetric(H))||(~isreal(H))`.
  - Representative operation: `(~issymmetric(H))||(~isreal(H))`.

## Parameters / inputs

- H -Hessian matrix to be regularised
- g -gradient computed at the same point as H
- data -diagnostic data structure

## Outputs

- H -regularised Hessian
- data -updated diagnostic data structure with
- data.count.rfo incremented by the number
- of RFO iterations taken

## Implementation structure

- RFO regularisation for Newton-Raphson Hessian and gradient
- pairs. Syntax:
- [H,data]=hessreg(spin_system,H,g,data)
- H -Hessian matrix to be regularised
- g -gradient computed at the same point as H
- data -diagnostic data structure
- H -regularised Hessian
- data -updated diagnostic data structure with
- data.count.rfo incremented by the number
- of RFO iterations taken
- Check consistency
- Set shorthands

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `logical()`, `cond()`, `speye()`, `issymmetric()`, `iscolumn()`.
