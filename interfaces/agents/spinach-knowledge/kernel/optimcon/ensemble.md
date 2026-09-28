# kernel/optimcon/ensemble.m

- Signature: `[traj_data,fidelity,gradient,hessian]=ensemble(waveform,spin_system)`

## Purpose

A parallel wrapper around GRAPE for ensemble optimal control optimisation. It handles multiple control power levels, resonance offsets, multistate transfers, and ensembles of drift Liouvillians.

## Physical / mathematical content

GRAPE propagates fidelity derivatives through a piecewise-constant pulse sequence for gradient-based optimisation of waveform samples.

## Numerical / algorithmic content

Parallel execution distributes ensemble cases across workers. The per-case propagation and derivative calculations run in `ens_block.m`; each worker sums its gradients, Hessians, and trajectories before the results are collected. This matters for large Spinach operators arising from basis expansion or powder or spatial lifting.

## Parameters / inputs

- `waveform` — control coefficients for each control operator, in rad/s.
- `spin_system` — spin system containing the ensemble control problem prepared by `optimcon.m`.

## Outputs

- `traj_data` — trajectory data for diagnostic plotting. Trajectories are returned in catalog order, or as an ensemble average when the `average` trajectory option is selected.
- `fidelity` — ensemble-averaged figure of merit for overlap between the current and desired state(s). With penalty methods, an array separates penalties from simulation fidelity.
- `gradient` — ensemble-averaged fidelity gradient with respect to the control sequence. With penalty methods, an array separates penalty gradients from the fidelity gradient.
- `hessian` — ensemble-averaged fidelity Hessian with respect to the control sequence. With penalty methods, an array separates penalty Hessians from the fidelity Hessian.

## Ensemble worker data flow

Cases in `spin_system.control.catalog` are assigned contiguous per-worker blocks by `optimcon.m` in `spin_system.control.worker_cases`. Each worker holds the common frozen problem and its block's drift generators as pool constants. At each objective evaluation, only the waveform and live control fields travel to the workers; per-case physics runs in `ens_block.m`. Call `ensemble` from the client using the same pool that was open when `optimcon.m` ran, because each worker holds only its own block.

[Source documentation](https://spindynamics.org/wiki/index.php?title=ensemble.m)