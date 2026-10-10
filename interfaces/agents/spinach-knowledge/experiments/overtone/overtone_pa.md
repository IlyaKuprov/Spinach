# experiments/overtone/overtone_pa.m

Source: [canonical MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/overtone/overtone_pa.m) · [Spinach wiki](https://spindynamics.org/wiki/index.php?title=overtone_pa.m)

Signature: `spectrum=overtone_pa(spin_system,parameters,H,R,K)`.

## Behaviour

This is an overtone soft-pulse/acquire wrapper. It obtains the source-defined reference frequency, lifts the supplied quadrupolar-channel `Lx` operator across `spc_dim`, applies one pulse to `rho0`, then calls `overtone_a` with the original `H`, `R`, and `K`. The pulse evolution combines the supplied Hamiltonian, relaxation, and kinetics as `H+1i*R+1i*K`.

The `average` branch forms an average pulse operator at `2*pi*(ovt_frq-rf_frq)` and applies its propagator for `rf_dur`. The `fplanck` branch calls `shaped_pulse_af` with the same combined evolution generator, `Lx`, the starting state, the offset in Hz, pulse power, and duration. The file has no separate MAS-rate or sample-orientation field; `spc_dim` is the documented Fokker–Planck spatial dimension.

## Inputs and limits

Required fields: `sweep` (two-element real frequency extent in Hz around the overtone frequency), `npoints` (positive integer), `spins` (one-element cell array), `spc_dim` (positive integer), `rho0`, `coil`, `Lx` (quadrupolar-nucleus X Zeeman operator), `rf_frq` (offset in Hz), `rf_pwr` (rad/s), `rf_dur` (positive scalar in seconds), and `method` (`'average'` or `'fplanck'`). The consistency check requires `H`, `R`, and `K` to be numeric matrices of matching dimensions.

The source comment states that relaxation must be present for the matrix inversion in `overtone_a` to converge and that `R` must not be thermalised. This is a documented requirement, not evidence of a successful run. The wrapper returns the helper's `spectrum` without reshaping it; the helper interprets `sweep` and `npoints` for acquisition.
