# experiments/hyperpol/topdnp.m

- Signature: `contact_curve=topdnp(spin_system,parameters,H,R,K)`
- Canonical MATLAB source: https://github.com/IlyaKuprov/Spinach/blob/main/experiments/hyperpol/topdnp.m

## Purpose and pulse cycle

Implements the time-optimised pulsed DNP loop described in the cited paper. The function combines context-supplied matrices as `L=H+1i*R+1i*K`, forms an electron (X)-directed microwave pulse by adding `2*pi*irr_powers*Ex`, then repeats pulse followed by delay. It records the detection-state overlap before the first pulse and after every completed loop. The supported propagation branches are `zeeman-hilb`, `zeeman-liouv`, and `sphten-liouv`.

## Inputs and output

`H`, `R`, and `K` are context-supplied matrices. Required fields are `irr_powers` (non-negative microwave amplitude in Hz), `rho0` (initial state), `coil` (detection state), `pulse_dur` and `delay_dur` (seconds), and `nloops` (positive integer). The `2*pi` factor converts the Hz amplitude in the driven generator to angular frequency; pulse and delay lengths are in seconds.

`contact_curve` is a `1 x (nloops+1)` row: the initial `coil`–state overlap, followed by one detected value after each pulse–delay cycle. It is not an FID, MRI image, or k-space trajectory. The source has no gradient or spatial-encoding input.

The numeric examples directly encoded by the source are the x-axis microwave term and the initial-plus-one-sample-per-loop convention; neither the source nor baseline page supplies a numerical parameter set or calculated experiment result.

## References

- TOP DNP paper: https://doi.org/10.1126/sciadv.aav6909
- Spin Dynamics Wiki: https://spindynamics.org/wiki/index.php?title=topdnp.m
