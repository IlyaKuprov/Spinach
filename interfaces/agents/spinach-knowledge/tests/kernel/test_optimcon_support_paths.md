# tests/kernel/test_optimcon_support_paths.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_optimcon_support_paths.m`
- Signature: `result=test_optimcon_support_paths()`
- Total lines: 352

## Purpose

Tests small optimal-control support paths. Syntax: result=test_optimcon_support_paths()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- Finite-difference discretisation appears in the implementation, so numerical accuracy depends on stencil order, boundary handling, and the balance between resolution and conditioning.
- The file also defines local helper function(s): `local_spin_system()`, `local_data()`, `local_penalty_term()`, `local_finite_grad()`, `local_finite_hess()`, `local_penalty_grad()`, `local_objective()`, `local_total_objective()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 16-17: Announce the test target; implemented by `fprintf('TESTING: Optimal-control support helper paths\n')`.
- Lines 19-22: State the support-path target of the test; implemented by `result=new_test_result('optimcon/support_paths', 'Optimal-control support helper paths', 'penalty(), trapdiff(), objeval(), and line-search helpers must pass small deter…`.
- Lines 24-25: Make a minimal quiet Spinach object for low-level helper calls; implemented by `spin_system=local_spin_system()`.
- Lines 27-29: Check the no-penalty path; implemented by `waveform=[-2.0 -0.5 0.25 2.0 0.75 -1.25; 1.5 -1.5 0.10 -0.1 0.50 1.25]`.
- Lines 38-39: Check the norm-square penalty against its closed form; implemented by `[pen_ns,grad_ns,hess_ns]=penalty(waveform,'NS',-1,1)`.
- Lines 50-51: Check the spillout penalty against explicit clipping residuals; implemented by `[pen_sns,grad_sns,hess_sns]=penalty(waveform,'SNS',-1,1)`.
- Lines 67-68: Check the derivative norm-square gradient by finite differences; implemented by `[pen_dns,grad_dns,hess_dns]=penalty(waveform,'DNS',-10,10)`.
- Lines 78-82: Check the amplitude spillout path on Cartesian controls; implemented by `amp_waveform=[ 0.50 2.00 -0.25 0.20; 0.50 0.50 1.50 -0.10; -0.40 0.20 1.20 -1.60; 0.30 -0.10 0.90 -0.40]`.
- Lines 108-109: Check trapezium derivative matrices against centred finite differences; implemented by `S=pauli(2)`.
- Lines 127-128: Check objective evaluation for fidelity, gradient, and Hessian collection; implemented by `data=local_data([2 1])`.
- Lines 153-154: Check direct line-search condition helpers; implemented by `spin_system.control.freeze=[]`.
- Lines 175-176: Check bracketing accepts a short ascent step on a concave quadratic; implemented by `data_line=local_data([2 1])`.
- Lines 191-192: Check sectioning recovers the bracketed maximum of the same quadratic; implemented by `a.alpha=0`.

### Key state/data transformations

- Lines 20-22: computes `result` using `result=new_test_result('optimcon/support_paths', 'Optimal-control support helper paths', 'penalty(), trapdiff(), objeval(), and line-search helpers must pass small deter…`.
- Lines 25: computes `spin_system` using `spin_system=local_spin_system()`.
- Lines 28-29: computes `waveform` using `waveform=[-2.0 -0.5 0.25 2.0 0.75 -1.25; 1.5 -1.5 0.10 -0.1 0.50 1.25]`.
- Lines 30: computes `[pen_none,grad_none,hess_none]` using `[pen_none,grad_none,hess_none]=penalty(waveform,'none',-1,1)`.
- Lines 39: computes `[pen_ns,grad_ns,hess_ns]` using `[pen_ns,grad_ns,hess_ns]=penalty(waveform,'NS',-1,1)`.
- Lines 40: computes `pen_ns_ref` using `pen_ns_ref=sum(waveform(:).^2)/size(waveform,2)`.
- Lines 41: computes `grad_ns_ref` using `grad_ns_ref=2*waveform/size(waveform,2)`.
- Lines 42: computes `hess_ns_ref` using `hess_ns_ref=2*eye(numel(waveform))/size(waveform,2)`.
- Lines 51: computes `[pen_sns,grad_sns,hess_sns]` using `[pen_sns,grad_sns,hess_sns]=penalty(waveform,'SNS',-1,1)`.
- Lines 52: computes `spill_hi` using `spill_hi=max(waveform-1,0)`.
- Lines 53: computes `spill_lo` using `spill_lo=min(waveform+1,0)`.
- Lines 54: computes `pen_sns_ref` using `pen_sns_ref=sum(spill_hi(:).^2)+sum(spill_lo(:).^2)`.
- Lines 56: computes `grad_sns_ref` using `grad_sns_ref=2*(spill_hi+spill_lo)/size(waveform,2)`.
- Lines 57: computes `spill_mask` using `spill_mask=(abs(spill_hi)+abs(spill_lo)>0)`.
- Lines 58: computes `hess_sns_ref` using `hess_sns_ref=2*diag(spill_mask(:))`.
- Lines 68: computes `[pen_dns,grad_dns,hess_dns]` using `[pen_dns,grad_dns,hess_dns]=penalty(waveform,'DNS',-10,10)`.
- Lines 69-70: computes `grad_dns_ref` using `grad_dns_ref=local_finite_grad(@(x)local_penalty_term(x,size(waveform),'DNS',-10,10), waveform(:))`.
- Lines 79-82: computes `amp_waveform` using `amp_waveform=[ 0.50 2.00 -0.25 0.20; 0.50 0.50 1.50 -0.10; -0.40 0.20 1.20 -1.60; 0.30 -0.10 0.90 -0.40]`.

### Local helper functions

- Line 211: `local_spin_system()` — `function spin_system=local_spin_system()`. Build a minimal quiet Spinach object for matrix-level helpers
  - Representative operation: `spin_system.sys.output='hush'`.
  - Representative operation: `spin_system.sys.enable={}`.
- Line 224: `local_data()` — `function data=local_data(x_shape)`. Build the optimisation workspace fields used by objeval()
  - Representative operation: `data.x_shape=x_shape`.
  - Representative operation: `data.count.fx=0`.
- Line 235: `local_penalty_term()` — `function term=local_penalty_term(x,wf_size,type,fb,cb)`. Evaluate only the scalar penalty value for finite differences
  - Representative operation: `term=penalty(reshape(x,wf_size),type,fb,cb)`.
- Line 243: `local_finite_grad()` — `function grad=local_finite_grad(fun,x)`. Compute a centred finite-difference gradient
  - Representative operation: `step_size=1e-6`.
  - Representative operation: `grad=zeros(size(x))`.
- Line 259: `local_finite_hess()` — `function hess=local_finite_hess(fun,x)`. Compute a centred finite-difference Hessian
  - Representative operation: `step_size=1e-6`.
  - Representative operation: `grad=fun(x)`.
- Line 276: `local_penalty_grad()` — `function grad=local_penalty_grad(x,wf_size,type,fb,cb)`. Evaluate only the penalty gradient for finite differences
  - Representative operation: `[~,grad]=penalty(reshape(x,wf_size),type,fb,cb)`.
  - Representative operation: `grad=grad(:)`.
- Line 285: `local_objective()` — `function [traj_data,fidelity,grad,hess]=local_objective(x,~)`. Define a two-channel objective with a penalty-like second component
  - Representative operation: `traj_data.marker='local_objective'`.
  - Representative operation: `target=[1; -2]`.
- Line 308: `local_total_objective()` — `function value=local_total_objective(x)`. Return the scalar objective assembled by objeval()
  - Representative operation: `target=[1; -2]`.
  - Representative operation: `residual=x-target`.

## Outputs

- result -regression test result with explanatory messages
- The test covers penalty functions, trapezium-product derivatives,
- objective-function collection, and small line-search helper paths.

## Implementation structure

- Tests small optimal-control support paths. Syntax:
- result=test_optimcon_support_paths()
- result -regression test result with explanatory messages
- The test covers penalty functions, trapezium-product derivatives,
- objective-function collection, and small line-search helper paths.
- Announce the test target
- State the support-path target of the test
- Make a minimal quiet Spinach object for low-level helper calls
- Check the no-penalty path
- Check the norm-square penalty against its closed form
- Check the spillout penalty against explicit clipping residuals
- Check the derivative norm-square gradient by finite differences

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `penalty()`, `trapdiff()`, `objeval()`, `local_spin_system()`, `test_close()`, `waveform()`, `spill_hi()`, `spill_lo()`, `spill_mask()`, `local_finite_grad()`, `local_penalty_term()`, `test_true()`, `grad_dns()`, `isequal()`, `amp_waveform()`.
