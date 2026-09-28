# experiments/nqr/nqr_pa.m

- Signature: `spectrum=nqr_pa(spin_system,parameters,H,R,K)`

## Purpose

Nuclear quadrupole resonance soft pulse-acquire experiment with idealised, infinite-bandwidth acquisition.

## Parameters / inputs

- `spin_system` — spin system supplied by the simulation context.
- `parameters.sweep` — two-element vector specifying the spectrum window extents in Hz, in ascending order.
- `parameters.npoints` — number of points in the spectrum; a positive integer.
- `parameters.rho0` — initial state.
- `parameters.coil` — detection state.
- `parameters.Lx`, `parameters.Ly` — operators used in the RF Hamiltonian; numeric matrices of equal size.
- `parameters.rf_frq` — RF irradiation frequency in Hz.
- `parameters.rf_pwr` — multiplier in rad/s of `Lx*cos(ωt)+Ly*sin(ωt)` in the RF Hamiltonian.
- `parameters.rf_dur` — pulse duration in seconds; non-negative.
- `H` — Hamiltonian matrix received from the context function.
- `R` — relaxation superoperator received from the context function.
- `K` — kinetics superoperator received from the context function.

## Outputs

- `spectrum` — spectrum of the specified initial state detected on the specified coil state within the requested frequency interval.

Relaxation must be present in the system dynamics for the matrix inversion to converge. The relaxation matrix `R` should **not** be thermalised.

## Implementation

The function converts the simulation to Liouville space, converts the RF operators to commutation superoperators when starting in `zeeman-hilb`, applies a soft off-resonance pulse with `shaped_pulse_af`, and performs frequency-domain acquisition with `slowpass`.

Source: <https://spindynamics.org/wiki/index.php?title=nqr_pa.m>