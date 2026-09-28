# experiments/nmr_solids/redor.m

- Signature: `redor_curve=redor(spin_system,parameters,H,R,K)`

## Purpose

Simulates a rotational-echo double-resonance (REDOR) experiment with ideal hard `pi` pulses on the dephasing channel. Called from the `singlerot` context, it returns the full echo, dephased echo, and their difference at requested numbers of rotor cycles. See https://doi.org/10.1016/0022-2364(89)90280-1 and https://doi.org/10.1006/jmre.2000.2128.

## Parameters / inputs

- `spin_system` — spin system; requires `sphten-liouv` formalism.
- `parameters.spc_dim` — Fokker–Planck spatial dimension, received from the context function.
- `parameters.spins` — observed and dephasing spins, respectively; for example, `{'13C','15N'}` for `13C{15N}` REDOR.
- `parameters.ncycles` — row vector of non-negative integer rotor-cycle counts in the REDOR evolution time.
- `parameters.rate` — MAS rate in Hz; a non-zero real scalar.
- `parameters.rho0` — initial state vector, usually transverse magnetisation on the observed spin.
- `parameters.coil` — detection state vector, usually on the observed spin.
- `parameters.decouple` — nuclei to decouple analytically during REDOR evolution, as a cell array of isotope strings; defaults to `{}` and must not include the REDOR working spins.
- `parameters.pulse_phase` — phases of the dephasing-channel `pi` pulses in radians; the vector cycles through the pulse train and defaults to `[0 pi/2]`.
- `H` — Hamiltonian matrix, received from the context function.
- `R` — relaxation superoperator, received from the context function.
- `K` — kinetics superoperator, received from the context function.

## Outputs

- `redor_curve(1,:)` — full echo `S0`, without dephasing pulses.
- `redor_curve(2,:)` — dephased echo `S`, with dephasing-channel `pi` pulses every half rotor period.
- `redor_curve(3,:)` — REDOR difference, `S0-S`.

## Implementation

The function forms `L=H+1i*R+1i*K` and applies any requested analytical decoupling. It propagates the reference trajectory for one rotor period per cycle; the dephased trajectory propagates for two half-periods, each followed by a phased dephasing-channel `pi` pulse. Detection uses `parameters.coil'*rho`, including at zero cycles, and the result selects the counts in `parameters.ncycles`.

https://spindynamics.org/wiki/index.php?title=redor.m