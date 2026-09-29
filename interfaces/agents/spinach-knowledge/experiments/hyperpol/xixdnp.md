# experiments/hyperpol/xixdnp.m

- Signature: `contact_curve=xixdnp(spin_system,parameters,H,R,K)`
- Context: call from a powder context; `H`, `R`, and `K` come from the context.
- Canonical implementation: [`experiments/hyperpol/xixdnp.m`](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/hyperpol/xixdnp.m)

## What it computes

This routine follows the TPPM DNP family and its X-inverse-X (XiX) special case. Starting from `parameters.rho0`, it applies a repeated two-pulse electron-microwave contact block and records the detected overlap after the initial state and after every block. The spin dynamics include the supplied Hamiltonian, relaxation, and kinetics; the routine itself does not solve a steady state or return an image. It produces a simulated contact curve, not a measured polarisation result.

## Pulse block and units

The source builds `L=H+1i*R+1i*K`, electron controls `Ex` and `Ey`, an on-axis drive `L1=L+2*pi*irr_powers*Ex`, and a second drive `L2=L+2*pi*irr_powers*(Ex*cos(phase)+Ey*sin(phase))`. The amplitude `irr_powers` is a nutation frequency in Hz, converted to angular frequency with `2*pi`; `phase` is in radians. Each block applies the +X pulse and then the phase-set pulse, each lasting `pulse_dur` seconds. No inter-pulse or shot-spacing delay is specified by this function. `nloops` is a positive integer number of two-pulse blocks.

The implementation supports `zeeman-hilb`, `zeeman-liouv`, and `sphten-liouv` formalisms. In Hilbert space it evolves the density operator as `P*rho*P'`; in Liouville space it applies the two propagators to the state vector in sequence.

## Inputs and output

Required fields are `irr_powers` (Hz), numeric initial state `rho0`, detection state `coil`, positive pulse duration `pulse_dur` (seconds), positive integer `nloops`, and real scalar `phase` (radians). `H`, `R`, and `K` must be numeric matrices compatible with the selected formalism. The returned row vector has `nloops+1` values: element 1 is `hdot(coil,rho0)`; each subsequent element is the detection overlap after one more complete pulse pair.

## References

- Source paper: [DOI 10.1021/jacs.1c09900](https://doi.org/10.1021/jacs.1c09900).
- Spin Dynamics Wiki: [`xixdnp.m`](https://spindynamics.org/wiki/index.php?title=xixdnp.m).
- MATLAB implementation: [`experiments/hyperpol/xixdnp.m`](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/hyperpol/xixdnp.m).
