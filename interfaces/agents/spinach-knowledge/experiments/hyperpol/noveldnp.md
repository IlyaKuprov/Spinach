# experiments/hyperpol/noveldnp.m

- Signature: `contact_curve=noveldnp(spin_system,parameters,H,R,K)`
- Canonical MATLAB source: https://github.com/IlyaKuprov/Spinach/blob/main/experiments/hyperpol/noveldnp.m

## Purpose and sequence

Computes a DNP contact-time observable for the pulsed solid effect or NOVEL. The supplied context matrices are combined as `L=H+1i*R+1i*K`; `K` is the kinetics superoperator, so it is not an MRI readout or a substitute for the spin Hamiltonian. The source builds electron-spin operators from the electron raising operator.

Both branches use a microwave spin-lock contact period along (-Y). With `parameters.flippulse=0`, propagation starts directly from `parameters.rho0` (the no-prepulse solid-effect branch). With `parameters.flippulse=1`, an electron (X)-axis 90-degree pulse is first applied for `parameters.pulse_dur`, then the (-Y) contact period is sampled (NOVEL). The microwave amplitude is given in Hz and multiplied by `2*pi` in the generator.

## Inputs and units

`H`, `R`, and `K` are context-supplied matrices. Required fields are `parameters.irr_powers` (non-negative microwave amplitude, Hz), `rho0` (initial state), `coil` (detection state), `timestep` (positive seconds), `nsteps` (positive integer), and `flippulse` (0 or 1). When `flippulse=1`, `pulse_dur` is also required and is a positive duration in seconds.

## Output and scope

`contact_curve` is the single-coil observable trace from the `evolution(...,'observable')` path: the initial value followed by the values at each of `nsteps` time steps. It is an observable curve, not a polarisation measurement, FID, image, or k-space array. This function has no gradient or spatial-encoding input.

Source-defined numeric choices include the 0/1 pulse switch and the 90-degree preparation followed by a 270-degree ((-Y)) spin-lock axis. The source and baseline page do not give an example parameter set or a computed numeric result.

## References

- NOVEL/solid-effect references retained from the source: https://doi.org/10.1016/0022-2364(88)90190-4 and https://doi.org/10.1063/1.5000528
- Spin Dynamics Wiki: https://spindynamics.org/wiki/index.php?title=noveldnp.m
