# experiments/nmr_liquids/hmbc.m

- Signature: `fid=hmbc(spin_system,parameters,H,R,K)`

## Purpose

Magnitude-mode HMBC pulse sequence.

## Implementation

The sequence applies an F2 90-degree pulse, evolves for `abs(1/(2*parameters.J))`, applies an F1 90-degree pulse, and evolves for `parameters.delta_b`. A difference of positive- and negative-angle F1 pulses selects coherence, followed by `+1` F1 coherence selection. F1 evolution is split around an F2 180-degree decoupling pulse; a final F1 pulse precedes F2 detection. The source cites 60 ms as the authors' recommended value for `parameters.delta_b`; the parameter is required and has no default in the function.

## Parameters / inputs

- `parameters.sweep` — `[F1 F2]` sweep widths in each frequency direction, Hz.
- `parameters.npoints` — `[F1 F2]` numbers of points in each time direction.
- `parameters.spins` — `{F1 F2}` nuclei, e.g. `{'15N','1H'}`.
- `parameters.J` — primary scalar coupling, Hz.
- `parameters.delta_b` — `delta_2` delay from the cited paper; the authors recommend `60e-3` seconds.
- `H` — Hamiltonian matrix, received from context function.
- `R` — relaxation superoperator, received from context function.
- `K` — kinetics superoperator, received from context function.

The implementation requires two positive sweep widths, two positive integer point counts, distinct spin isotopes present in the system, non-zero real `J`, and positive real `parameters.delta_b`.

## Outputs

- `fid` — free induction decay for magnitude-mode processing.

Natural-abundance experiments should make use of the isotope-dilution functionality; see `dilute.m`.

## References

- [HMBC reference](https://doi.org/10.1021/ja00268a061)
- [HMBC reference](https://doi.org/10.1016/0022-2364(88)90172-2)
- [Spin Dynamics Wiki: `hmbc.m`](https://spindynamics.org/wiki/index.php?title=hmbc.m)
