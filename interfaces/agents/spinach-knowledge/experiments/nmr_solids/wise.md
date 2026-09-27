# experiments/nmr_solids/wise.m

- Signature: `fid=wise(spin_system,parameters,H,R,K)`

## Purpose

WISE (WIdeline SEparation) is a powder MAS heteronuclear correlation experiment. In the common 1H–13C implementation, molecular dynamics information is contained in 1H line shapes separated in the second dimension by 13C chemical shifts. See https://doi.org/10.1021/ma00038a037.

## Parameters / inputs

- `spin_system` — spin system; this implementation requires the `sphten-liouv` formalism.
- `parameters.spins` — two working spin isotopes, e.g. `{'1H','13C'}`.
- `parameters.hi_pwr` — amplitude of high-power pulses on the high-gamma channel, Hz; a positive real scalar.
- `parameters.cp_pwr` — RF amplitudes on the two channels during CP contact, Hz; a two-element row vector of positive values.
- `parameters.cp_dur` — CP contact time duration, s; a positive real scalar.
- `parameters.rho0` — initial state.
- `parameters.coil` — detection state.
- `parameters.sweep` — sweep widths for F1 and F2, Hz; a two-element row vector of positive values.
- `parameters.npoints` — numbers of points in F1 and F2; a two-element row vector of positive integers.
- `H` — Hamiltonian matrix, received from the context function.
- `R` — relaxation superoperator, received from the context function.
- `K` — kinetics superoperator, received from the context function. `H`, `R`, and `K` must have the same dimensions.

## Outputs

- `fid.sin`, `fid.cos` — sine and cosine components of the States quadrature.

## Implementation summary

The sequence pre-saturates the second spin channel, applies separate high-power 90-degree pulses to the first channel for the cosine and sine pathways, and evolves both pathways during F1. Cross-polarization follows, with RF terms on both channels for `parameters.cp_dur`. The first channel is then wiped and decoupled before F2 acquisition through `parameters.coil`.

Source reference: <https://spindynamics.org/wiki/index.php?title=wise.m>