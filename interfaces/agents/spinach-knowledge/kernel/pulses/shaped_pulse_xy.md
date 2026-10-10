# kernel/pulses/shaped_pulse_xy.m

- MATLAB source: [kernel/pulses/shaped_pulse_xy.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/pulses/shaped_pulse_xy.m)
- Spinach wiki: [shaped_pulse_xy.m](https://spindynamics.org/wiki/index.php?title=shaped_pulse_xy.m)
- Signature: `[rho,traj,P]=shaped_pulse_xy(spin_system,drift,controls,amplitudes,slice_durs,rho,method)`

## Purpose

Applies a Cartesian shaped pulse on supplied control operators while the drift Liouvillian continues to act. The drift should include the transmitter offset, if present; controls may include spatial operators such as gradients and diffusion.

## Slice construction

For piecewise-constant quadrature, each slice uses `drift + sum(amplitudes{k}(n)*controls{k})`. For piecewise-linear quadrature, the source uses the left- and right-edge amplitude values for each control and combines the edge generators with `isergen`, its two-point, second-order Lie quadrature. Slice durations are in seconds; control-amplitude values are in radians per second. Each amplitude vector has one entry per slice for piecewise-constant quadrature and one extra entry for piecewise-linear quadrature.

The propagation options are `'expv-pwc'`, `'expv-pwl'`, `'expm-pwc'`, `'expm-pwl'`, `'evol-pwc'`, and `'evol-pwl'`. Here `expv` is the Krylov method, `expm` uses explicit matrix exponentiation, and `evol` calls Spinach evolution; the source advises against `evol` unless there is a specific reason. The `PWL` suffix selects the two-point quadrature; `PWC` uses a constant slice operator.

## Outputs and formalism

- `rho` — final state vector or stack of states.
- Optional `traj` — a `1 x (nsteps+1)` cell array whose first element is the initial condition.
- Optional `P` — effective pulse propagator; the source describes its construction as expensive. In Hilbert-space density-matrix form, apply it on both sides: `P*rho_initial*P'`.

The source applies each slice propagator by left multiplication for `sphten-liouv`, `zeeman-liouv`, and `zeeman-wavef` formalisms; for `zeeman-hilb` it applies `P*rho*P'`.
