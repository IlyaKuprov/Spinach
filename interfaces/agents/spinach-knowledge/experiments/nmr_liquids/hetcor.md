# experiments/nmr_liquids/hetcor.m

- Signature: `fid=hetcor(spin_system,parameters,H,R,K)`

## Purpose

Magnitude-mode HETCOR pulse sequence. The fixed transfer delays are `delta_2=1/(2J)` and `delta_3=1/(3J)`, using the absolute value of `J`.

## Implementation

The initial F1 pulse selects coherence, followed by F1 evolution with an F2 refocusing pulse. After `delta_2`, simultaneous F1/F2 pulses select the transfer pathway; `delta_3` follows. The requested F1 channels are decoupled before F2 detection, producing a two-dimensional free induction decay for magnitude-mode processing. The code sets `delta_2=abs(1/(2*parameters.J))` and `delta_3=abs(1/(3*parameters.J))`.

## Parameters / inputs

- `parameters.sweep` — `[F1 F2]` sweep widths in the two frequency directions, Hz.
- `parameters.npoints` — `[F1 F2]` numbers of points in the two time directions in `fid`.
- `parameters.spins` — `{F1 F2}` nuclei (e.g. `{'1H','13C'}`).
- `parameters.decouple` — list of nuclei that detection-time decoupling should be applied to, a cell array of strings (e.g. `{'1H','15N'}`).
- `parameters.J` — working scalar coupling, Hz.
- `H` — Hamiltonian matrix, received from context function.
- `R` — relaxation superoperator, received from context function.
- `K` — kinetics superoperator, received from context function.

The implementation requires two positive sweep widths, two positive integer point counts, two different isotopes present in the system, a cell array of present decoupling isotopes, and non-zero real `J`.

## Outputs

- `fid` — two-dimensional free induction decay for magnitude-mode processing.

Natural-abundance experiments should use the isotope-dilution functionality; see `dilute.m`.

## References

- [HETCOR reference](https://doi.org/10.1016/0022-2364(81)90272-9)
- [Spin Dynamics Wiki: `hetcor.m`](https://spindynamics.org/wiki/index.php?title=hetcor.m)
