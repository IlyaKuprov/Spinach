# experiments/hp_acquire.m

- Signature: fid=hp_acquire(spin_system,parameters,H,R,K)

## Purpose and physical scope

Runs a user-specified hard-pulse/acquisition sequence. The source describes a standard pulse-acquire experiment; the function name does not establish that the initial state was produced by DNP or is hyperpolarised. Hyperfine couplings and relaxation affect the calculated signal only through the supplied system and context matrices. This is not an ESEEM or ENDOR experiment by itself.

## Inputs and parameters

H, R and K are context-supplied Hamiltonian, relaxation and kinetics matrices with matching dimensions. The routine combines them as H+1i*R+1i*K.

- parameters.sweep: acquisition sweep width in Hz.
- parameters.npoints: positive integer number of FID samples.
- parameters.rho0: initial state.
- parameters.coil: detection state.
- parameters.pulse_op: user-supplied pulse operator.
- parameters.pulse_angle: hard-pulse rotation angle in radians.
- parameters.decouple: required cell array of spin labels to decouple, or an empty cell array. The source example is {'15N','13C'}; decoupling is documented for sphten-liouv formalism.
- parameters.echo_time: optional evolution interval for echo detection. The source comment does not state its unit.
- parameters.echo_oper and parameters.echo_angle: echo pulse operator and angle, used when echo_time is supplied; angle is in radians.

## Sequence and output

The function applies pulse_op to rho0 with step and pulse_angle. If echo_time is present, it evolves for that interval, applies echo_oper through a second step, and evolves for the same interval again. It then applies the requested decoupling and calls evolution in observable mode with parameters.coil. The acquisition interval is 1/sweep, with npoints-1 evolution intervals, so fid contains the acquired observable samples; only fid is returned, not a separate time axis.

## Limits and interpretation

This routine provides pulse, optional echo, decoupling, and FID propagation under the supplied H, R and K. It does not create a polarisation mechanism, define a measured polarisation level, or report a separate echo amplitude. The source gives no numerical pulse angle, sweep width, or echo-time example; the isotope labels above are the source's actual decoupling example.

Source: https://github.com/IlyaKuprov/Spinach/blob/main/experiments/hp_acquire.m
