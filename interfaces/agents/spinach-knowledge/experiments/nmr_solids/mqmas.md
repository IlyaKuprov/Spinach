# experiments/nmr_solids/mqmas.m

Source: [canonical MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/nmr_solids/mqmas.m) · [Spin Dynamics Wiki](https://spindynamics.org/wiki/index.php?title=mqmas.m)

## Purpose and inputs

A rotor-synchronous MQMAS pulse sequence returning a 2D amplitude-mode free-induction decay. Call it through the `singlerot.m` context, which supplies `H`, `R`, and `K`. Signature: `fid=mqmas(spin_system,parameters,H,R,K)`. The three superoperators must be square numeric matrices of equal size.

- `parameters.spins` is a one-element cell containing an isotope string present in the system; `spc_dim` is a positive integer spatial dimension.
- `pulse_dur` is a two-element vector of non-negative durations in seconds. `pulse_amp` is a two-element real vector in rad/s; the code adds these values directly to the pulse generators.
- `mq_order` is an integer coherence order. `rho0` and `coil` are required initial and detection states.
- `npoints` is a two-element vector of positive integers. `rate` is a non-zero real MAS rate in Hz, and `sweep` must be positive and equal to `abs(rate)`. Both dimensions are sampled stroboscopically at `1/abs(rate)`.
- `decouple` is a cell array of isotope strings, possibly empty; listed isotopes must occur in the spin system. Analytical decoupling is restricted to `sphten-liouv` when the list is non-empty.

## Sequence outline

The source forms `L=H+1i*R+1i*K` and applies the decoupling configuration. It constructs `Lx` from the active spin’s `L+` operator, applies pulse 1, and selects `mq_order` coherence. The indirect trajectory advances for `npoints(1)-1` rotor periods. Pulse 2 is followed by selection of +1 coherence; direct acquisition then uses `coil` for `npoints(2)-1` intervals. The returned `fid` contains the two sampled dimensions in amplitude mode.

This mapping is MQMAS; the source does not describe an overtone cross-polarisation transfer block.