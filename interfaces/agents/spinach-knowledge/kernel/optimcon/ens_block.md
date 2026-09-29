# kernel/optimcon/ens_block.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/optimcon/ens_block.m) · [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=ens_block.m)

- Signature: `[traj,fid,grad,hess]=ens_block(spin_system,drifts,control,block,waveform,n_outputs)`

## Role and case mapping

This is the worker-side evaluator called by `ensemble.m`. For each case assigned to `block`, catalog columns select the state-target pair, drift generator, power level, offset combination, phase-cycle row, and distortion row. The selected initial and target states, drift and waveform are passed to the GRAPE engine chosen for the spin-system formalism.

Offsets are supplied in Hz and enter the drift as angular-frequency terms `2*pi*offset*off_op`. A phase-cycle row phases the initial and target states and rotates the control channels. The selected power level scales the rotated waveform; the case's distortion functions are then applied in column order. When derivatives are requested, their Jacobians are composed and used to pull the GRAPE gradient back to the input waveform. Frozen waveform samples have zero gradient entries.

## Outputs and limits

- `traj`: a cell array of trajectory structures, one per case in block order (empty for an empty block). With `'average'` in `control.traj_opts`, a nonempty block is collapsed to one structure whose `forward` member is the sum of that block's forward trajectories; `ensemble.m` forms the catalog-wide mean.
- `fid`: one fidelity per case, `1 x n_block`.
- `grad`: block-summed case gradients, a `(ncontrols*nsteps) x 1` column when requested.
- `hess`: block-summed case Hessians, flattened as a `((ncontrols*nsteps)^2) x 1` column when requested.

The waveform is real numeric, in rad/s, with `ncontrols` rows and `pulse_nsteps` columns for the rectangle integrator or one extra column for the trapezium integrator. The implementation checks type and shape but does not explicitly test sample finiteness. Hessians are supported only with the rectangle integrator and no waveform distortions. The Hessian pullback accounts for phase cycling and power scaling; frozen input coordinates have zero rows and columns.

## Inputs

`spin_system` and `drifts` are worker-resident frozen problem data published by `optimcon.m`; `control` is the live control structure supplied by `ensemble.m`; `block` is a valid positive integer worker-block index. `n_outputs` requests fidelity only, gradient, or Hessian through `ensemble.m` (2, 3, or 4 outputs respectively).
