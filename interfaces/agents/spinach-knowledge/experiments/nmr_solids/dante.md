# experiments/nmr_solids/dante.m

Source: [canonical MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/nmr_solids/dante.m) · [Spin Dynamics Wiki](https://spindynamics.org/wiki/index.php?title=dante.m)

## Purpose

A rotor-synchronous DANTE pulse train, normally called from the `singlerot.m` context that supplies the Hamiltonian, relaxation, and kinetics superoperators. This source describes DANTE; it should not be conflated with a REDOR sequence.

## Inputs

Signature: `fid=dante(spin_system,parameters,H,R,K)`.

- `H`, `R`, and `K` must be numeric matrices of equal size.
- `parameters.spins` is a one-element cell array naming an isotope present in the spin system. `parameters.decouple` is a cell array of isotope strings to decouple; it may be empty, and listed isotopes must be present.
- `parameters.pulse_dur` is a positive pulse duration in seconds; `pulse_amp` is a real RF amplitude in Hz; `pulse_num` and `n_periods` are positive integers.
- `parameters.rate` is a non-zero real rotor rate in Hz. `parameters.sweep` is a positive real acquisition sweep width in Hz, and `npoints` is a positive integer.
- `parameters.spc_dim` is a positive integer Fokker–Planck spatial dimension; `rho0` and `coil` are required initial and detection states.

## Propagation and output

The implementation forms `L=H+1i*R+1i*K`, builds `Lx` from the `L+` operator of `parameters.spins{1}`, and applies the requested decoupling. Each rotor period has duration `abs(1/rate)` and is divided into `pulse_num` slots. In every slot, the working spin evolves for `pulse_dur` under `L+2*pi*pulse_amp*Lx`, then freely under `L` for the rest of that slot. The pulse must fit its slot. After `n_periods`, the function acquires `fid` using `coil`, dwell `1/sweep`, and `npoints-1` propagation steps.

The output is the free-induction decay `fid`. This source provides no REDOR-specific echo/recoupling step or protein-specific transfer block.