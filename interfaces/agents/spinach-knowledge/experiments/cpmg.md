# experiments/cpmg.m

- MATLAB implementation: [experiments/cpmg.m](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/cpmg.m)

- Signature: `fid=cpmg(spin_system,parameters,H,R,K)`

## Purpose and signal

This routine generates a CPMG echo train from a caller-supplied initial state and records its projection onto a supplied detection state. It applies a nominal `pi/2` excitation, samples the first half-echo, then repeats the refocusing `pi` pulse and full-echo acquisition for `parameters.nloops` loops. The implementation does not describe a vendor-specific Bruker acquisition, receiver dead time, or phase cycle; those are not part of this function's contract.

## Inputs and pulse preparation

`H`, `R`, and `K` are numeric matrices of equal size and are combined as `H + 1i*R + 1i*K`. The parameter structure supplies `rho0` (initial state), `coil` (detection state), `pulse_op` (pulse operator), `nloops` (positive integer loop count), `timestep` (positive scalar propagation step), and `npoints` (number of steps in a half-echo). The routine lifts `pulse_op` with `kron(speye(parameters.spc_dim),parameters.pulse_op)`; `parameters.spc_dim` must therefore be present and consistent with the operator/state representation. The source does not assign a unit to `timestep`.

## Output

`fid` is a row vector of coil-detected samples: the initial half-echo followed by the sampled echoes. The detection is the conjugate-transpose projection `parameters.coil' * trajectory`, so complex signal values are retained. `npoints` sets the initial half-echo sampling; each repeated echo is propagated over `2*parameters.npoints-1` steps.

[Source page](https://spindynamics.org/wiki/index.php?title=cpmg.m)
