# experiments/nmr_liquids/hmqc.m

- Signature: `fid=hmqc(spin_system,parameters,H,R,K)`

## Purpose

Magnitude-mode HMQC pulse sequence.

## Implementation

The kernel starts from F2 `Lz` magnetisation, applies an F2 90-degree pulse and evolves for `abs(1/(2*parameters.J))`. An F1 90-degree pulse and `+1` F1 coherence selection precede indirect evolution. The F1 evolution is split around midpoint 180-degree refocusing pulses on `parameters.decouple_f1`; an F1 pulse and the second fixed coupling delay follow. Channels in `parameters.decouple_f2` are decoupled before F2 detection. The resulting free induction decay is for magnitude-mode processing.

## Parameters / inputs

- `parameters.sweep` — `[F1 F2]` sweep widths in each frequency direction, Hz.
- `parameters.npoints` — `[F1 F2]` numbers of points in each time direction.
- `parameters.spins` — `{F1 F2}` nuclei, e.g. `{'15N','1H'}`.
- `parameters.decouple_f2` — nuclei to decouple in F2, e.g. `{'15N','13C'}`.
- `parameters.decouple_f1` — nuclei that receive midpoint 180-degree refocusing pulses in F1, e.g. `{'1H'}`.
- `parameters.J` — primary scalar coupling, Hz.
- `H` — Hamiltonian matrix, received from context function.
- `R` — relaxation superoperator, received from context function.
- `K` — kinetics superoperator, received from context function.

The implementation requires two positive sweep widths, two positive integer point counts, distinct spin isotopes present in the system, valid decoupling isotope lists, and non-zero real `J`.

## Outputs

- `fid` — free induction decay for magnitude-mode processing.

Natural-abundance experiments should make use of the isotope-dilution functionality; see `dilute.m`.

## References

- [HMQC reference](https://doi.org/10.1016/0022-2364(83)90241-X)
- [Spin Dynamics Wiki: `hmqc.m`](https://spindynamics.org/wiki/index.php?title=hmqc.m)
