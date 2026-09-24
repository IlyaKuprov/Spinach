# kernel/optimcon/ens_block.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/optimcon/ens_block.m`
- Signature: `[traj,fid,grad,hess]=ens_block(spin_system,drifts,control,block,waveform,n_outputs)`
- Total lines: 208

## Purpose

Fidelity, gradient, and Hessian contributions of one block of ensemble cases, evaluated on the parallel pool worker that holds the drift generators of that block. This function is called by ensemble.m inside its spmd block; the per-case physics (phase cycle, offsets, power level, waveform distortions, GRAPE) is applied here. Syntax: `[traj,fid,grad,hess]=ens_block(spin_system,drifts,control,block,waveform,n_outputs)`.

## Physical / mathematical content

Each ensemble member has its own initial/target state pair, drift generator, power level, resonance offsets, phase-cycle step, and distortion model. Offsets specified in Hz contribute angular frequencies through the factor `2*pi`. Phase cycles rotate paired control channels and phase the states; power scaling and distortions determine the actual field experienced by the spins. The chosen GRAPE engine evaluates that member in the configured Hilbert or Liouville representation.

## Numerical / algorithmic content

For a Cartesian input waveform, the gradient is the pullback of the physical-field derivative through the distortion Jacobian, phase rotation, and power scaling. Only then are frozen Cartesian input entries zeroed. Supported Hessians undergo the corresponding two-sided phase/power transformation and have frozen rows and columns zeroed; distortion Hessians remain unavailable. An outer coordinate wrapper must defer its own mask: `grape_curv` passes an empty Cartesian mask and freezes only after its curvilinear pullback. Direct GRAPE engine calls retain their separate existing behaviour.

Worker-local drift storage avoids duplicating the complete drift ensemble on every worker. Fidelities remain per-case, while gradients and Hessians are block sums for reduction and ensemble averaging by `ensemble`. With `traj_opts` containing `average`, one summed forward trajectory is returned per non-empty block; otherwise trajectories retain block order. These sums also support Hilbert-space cell trajectories. An empty block contributes no trajectories.

## Syntax

`[traj,fid,grad,hess]=ens_block(spin_system,drifts,control,block,waveform,n_outputs)`

## Parameters / inputs

- `spin_system`: frozen problem published by `optimcon`, with the drift generators removed.
- `drifts`: cell array of drift generators populated at the indices used by this block's cases.
- `control`: live client-side control structure.
- `block`: index into `spin_system.control.worker_cases`.
- `waveform`: coefficients for the control operators, with controls in rows and time samples in columns; power scaling converts them to the fields used by the engines.
- `n_outputs`: requested output count: two for fidelity, three for the gradient, and four for the Hessian.

## Outputs

- `traj`: trajectories in block order, or one forward-trajectory sum when `traj_opts` contains `average`; an empty block returns an empty cell.
- `fid`: a `1 x n_block` array of per-case fidelities.
- `grad`: summed case gradients as a `ncontrols*nsteps x 1` column, empty unless `n_outputs>2`.
- `hess`: summed case Hessians as a `(ncontrols*nsteps)^2 x 1` column, empty unless `n_outputs>3`.

## Header notes

The worker receives the frozen spin-system description and its local drift generators together with live control data, the block identifier, waveform, and requested output count.
