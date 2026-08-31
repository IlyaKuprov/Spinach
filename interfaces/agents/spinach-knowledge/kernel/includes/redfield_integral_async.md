# kernel/includes/redfield_integral_async.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/includes/redfield_integral_async.m`
- Signature: `job_id=brw_compute_kernel(spin_system,w,job_id,upper_lim)`
- Total lines: 193

## Purpose

Bloch-Wangsness-Redfield and Nakajima-Zwanzig integral evaluation, the asynchronous parallel path. This include is called from within the relaxation.m theory blocks and follows the notation used in IK's paper: with the difference that the numerical quadrature method proposed there has been superceded by the much faster auxiliary matrix me- thod described in: The calling theory block must set rlx_onshell (true selects

## Physical / mathematical content

- Include scripts and shared setup fragments. These files implement tightly scoped runtime setup, parallel profiling, resource guards, or shared kernels included by other Spinach routines.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- Orientation or trajectory averaging is performed numerically, so grid design, weights, and integration error control matter directly to accuracy and runtime.

## Numerical / algorithmic content

- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 165-166: Get data from pool ValueStore; implemented by `store=getCurrentValueStore()`.
- Lines 173-174: Compute the integral and do an intermediate clean-up; implemented by `R=-w*ABCD{1}*expmint(spin_system,ABCD{2},ABCD{3},ABCD{4},upper_lim)`.
- Lines 177-178: Convert CSR sparse format into XYZ format; implemented by `[row,col,val]=find(R); clear('spin_system','R')`.
- Lines 180-182: Send to the pool ValueStore and delete local copy; implemented by `put(store,{['redfield_int_batch_' num2str(job_id)]}, {[row col val]})`.

### Key state/data transformations

- Lines 166: computes `store` using `store=getCurrentValueStore()`.
- Lines 167-170: computes `store_keys` using `store_keys={['brw_integrator_batch_' num2str(job_id) '_A'], ['brw_integrator_batch_' num2str(job_id) '_B'], ['brw_integrator_batch_' num2str(job_id) '_C'], ['brw_integra…`.
- Lines 171: computes `ABCD` using `ABCD=get(store,store_keys); remove(store,store_keys)`.
- Lines 174: computes `R` using `R=-w*ABCD{1}*expmint(spin_system,ABCD{2},ABCD{3},ABCD{4},upper_lim)`.
- Lines 175: computes `clear('ABCD'); R` using `clear('ABCD'); R=clean_up(spin_system,R,1e-2*spin_system.tols.rlx_zero)`.
- Lines 178: computes `[row,col,val]` using `[row,col,val]=find(R); clear('spin_system','R')`.

## Implementation structure

- Bloch-Wangsness-Redfield and Nakajima-Zwanzig integral evaluation,
- the asynchronous parallel path. This include is called from within
- the relaxation.m theory blocks and follows the notation used in
- IK's paper:
- with the difference that the numerical quadrature method proposed
- there has been superceded by the much faster auxiliary matrix me-
- thod described in:
- The calling theory block must set rlx_onshell (true selects the
- back-rotated kernel that reduces to Redfield theory at zero shift,
- false the resolvent kernel of Nakajima-Zwanzig theory) and
- rlx_shift (the Laplace evaluation point, Hz); Redfield theory is
- the on-shell form at zero shift.

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `cheap_norm()`, `corrfun()`, `report()`, `num2str()`, `clean_up()`, `speye()`, `gcp()`, `parfeval()`, `clear()`, `exist()`, `fetchNext()`, `rethrow()`, `cell2mat()`, `get()`, `remove()`, `XYZ()`.
