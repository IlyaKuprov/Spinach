# experiments/bruker/hsqcetgp.m

- Signature: `fid=hsqcetgp(spin_system,parameters,H,R,K)`

## Purpose

Simulate the echo/antiecho, gradient-selected HSQC sequence. Gradient selection is represented analytically by coherence-order selection rather than explicit gradient pulses. The implementation follows the Bruker `hsqcetgp` pulse program and standard HSQC sequence.

## Physical / mathematical content

The sequence starts from longitudinal magnetisation on the F2 spin and uses two INEPT periods of `abs(1/(2*J))/2`, separated by simultaneous F1/F2 refocusing pulses. A proton trim pulse and transfer pulses precede F1 evolution in two halves with configured refocusing pulses. Opposite F1 coherence orders define the echo and antiecho pathways; refocusing and back-transfer then lead to F2 coherence selection, optional F2 decoupling, and detection.

## Numerical / algorithmic content

The code forms `L = H + 1i*R + 1i*K`, uses dwell times `1/sweep(1)` and `1/sweep(2)`, and delegates propagation, pulses, and coherence selection to Spinach's `evolution`, `step`, and `coherence` routines. The function requires the `sphten-liouv` formalism.

## Syntax

```matlab
fid=hsqcetgp(spin_system,parameters,H,R,K)
```

## Parameters / inputs

- `parameters.sweep`: two positive real sweep widths, `[F1 F2]`, Hz.
- `parameters.npoints`: two positive integer point counts, `[F1 F2]`.
- `parameters.spins`: two different isotope strings, `{F1 F2}`, present in the spin system (e.g. `{'13C','1H'}`).
- `parameters.decouple_f2`: cell array of isotopes to decouple in F2 (e.g. `{'15N','13C'}`).
- `parameters.decouple_f1`: cell array of isotopes receiving midpoint 180-degree refocusing pulses in F1 (e.g. `{'1H','15N'}`); it must not include the active F1 isotope.
- `parameters.J`: non-zero real scalar working scalar coupling, Hz.
- `parameters.trim_angle`: finite real proton trim-pulse angle, radians.
- `H`: Hamiltonian matrix.
- `R`: relaxation superoperator.
- `K`: kinetics superoperator. These matrices must have matching dimensions.

## Outputs

- `fid.pos` and `fid.neg`: echo and antiecho signal components.

For natural-abundance simulations, the source recommends isotope dilution; see [`dilute.m`](https://spindynamics.org/wiki/index.php?title=dilute.m).

## References and links

- [HSQC sequence reference](https://doi.org/10.1016/0009-2614(80)80041-8)
- [HSQC reference](https://doi.org/10.1002/cmr.a.10095)
- [Spinach documentation for `hsqcetgp.m`](https://spindynamics.org/wiki/index.php?title=hsqcetgp.m)
