# kernel/optimcon/ensemble.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/optimcon/ensemble.m) · [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=ensemble.m)

- Signature: `[traj_data,fidelity,gradient,hessian]=ensemble(waveform,spin_system)`

## Purpose and objective

Evaluates the GRAPE figure of merit for every selected case in the ensemble catalog, distributes cases across the workers, and reduces their outputs on the client. Each case is passed the configured `control.fidelity` method. The returned fidelity is the dot product of the raw fidelities in catalog order with the normalised case weights. Workers apply the same weights to the transformed derivatives before reduction. Without explicit drift weights, this is the arithmetic mean over the retained cases. The gradient has the waveform shape and the Hessian is a square matrix with one row and column per waveform sample.

The cases are defined by the catalog built in `optimcon.m`: state-target pairs, drift generators, power levels, resonance-offset combinations, phase-cycle rows and distortion functions. Per-case transformations and the chain rule for distortions are handled by `ens_block.m`.

## Inputs and constraints

- `waveform`: real numeric control coefficients in rad/s, with `ncontrols` rows. It has `pulse_nsteps` columns for the rectangle integrator and `pulse_nsteps + 1` for the trapezium integrator. The wrapper checks type and shape but does not explicitly test sample finiteness.
- `spin_system`: spin system whose ensemble problem has been prepared by `optimcon.m`. If `control.return_traj` is absent, the wrapper sets it to `false`.

Call on the same parallel pool that was open when `optimcon.m` ran. Inside a `parfor` body on a worker, where no pool is visible, `optimcon.m` must run in the same loop body; the ensemble cases are then evaluated serially on that worker. The function rejects a changed pool, modified frozen generators/operators, or a changed ensemble composition. A Hessian request requires the rectangle integrator and no waveform distortions.

## Outputs

- `traj_data`: an `n_cases x 1` cell array of per-case trajectory structures in catalog order. If `'average'` is in `control.traj_opts`, returns a one-element cell containing a structure whose `forward` member is the weighted ensemble-mean forward trajectory.
- `fidelity`: scalar mean figure of merit over catalog cases.
- `gradient`: mean fidelity derivative with respect to waveform samples, with the same dimensions as `waveform`, returned when requested.
- `hessian`: mean fidelity Hessian with dimensions `numel(waveform) x numel(waveform)`, returned when requested and supported.

For gate control, fidelity is the weighted mean of the squared normalised trace overlaps, not the square of a mean trace. The trajectory structures carry no density trajectory; no gate Hessian is returned.
