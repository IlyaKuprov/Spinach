# kernel/pulses/shaped_pulse_xy.m

- Signature: `[rho,traj,P]=shaped_pulse_xy(spin_system,drift,controls,amplitudes,slice_durs,rho,method)`

## Purpose

Apply a Cartesian shaped pulse on specified control operators while the drift generator continues to act, including any transmitter offset.

## Physical / mathematical content

A slice generator is the drift plus the amplitude-weighted control operators. Propagation uses `exp(-1i*H*dt)`: wavefunction and Liouville states evolve by left multiplication, while Hilbert density matrices evolve by `P*rho*P'`. Piecewise-linear pulses use the two-point Lie quadrature implemented by `isergen`.

## Numerical / algorithmic content

`expv-pwc` and `expv-pwl` use `step` for the state; these action-exponential methods are usually preferable for one or two outputs. `expm-pwc` and `expm-pwl` construct explicit slice propagators and are usually preferable for three outputs. `evol-pwc` and `evol-pwl` call `evolution`; choose them only for a specific reason.

For all six methods in `zeeman-hilb`, the third output is the one-sided ordered product of slice propagators. Requesting it with `expv-*` or `evol-*` explicitly constructs those propagators without changing the selected two-sided state-propagation method. Propagator accumulation uses the configured clean-up tolerance. GPU arithmetic is selected by `'gpu'` in `spin_system.sys.enable`; returned states, trajectories, and propagators are gathered to host memory.

## Syntax

```matlab
[rho,traj,P]=shaped_pulse_xy(spin_system,drift,controls,...
                            amplitudes,slice_durs,rho,method)
```

## Parameters / inputs

- `spin_system`: Spinach system carrying the formalism, numerical tolerances, and execution options.
- `drift`: background Hamiltonian or Liouvillian, including the transmitter offset when needed.
- `controls`: cell array of control operators, one per channel; operators may include spatial degrees of freedom such as gradients and diffusion.
- `amplitudes`: cell array of amplitude vectors in rad/s, one per control channel.
- `slice_durs`: pulse-slice durations in seconds. PWC amplitudes have one element per slice; PWL amplitudes have one extra element for the edge values.
- `rho`: initial state vector or bookshelf of vectors in one-sided formalisms; density matrix in `zeeman-hilb`.
- `method`: one of `expv-pwc`, `expv-pwl`, `expm-pwc`, `expm-pwl`, `evol-pwc`, or `evol-pwl`.

## Outputs

- `rho`: final state with the same physical interpretation as the input.
- `traj`: `1 x (nsteps+1)` cell array of states, including the initial condition.
- `P`: effective pulse propagator, expensive to construct and best avoided unless needed. In Hilbert space reuse it as `P*rho_initial*P'`, not `P*rho_initial`.

## Header notes

The source links to https://spindynamics.org/wiki/index.php?title=shaped_pulse_xy.m .
