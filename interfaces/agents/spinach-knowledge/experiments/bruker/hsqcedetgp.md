# experiments/bruker/hsqcedetgp.m

- Signature: `fid=hsqcedetgp(spin_system,parameters,H,R,K)`

## Purpose

Simulate the echo/antiecho, gradient-selected multiplicity-edited HSQC sequence. Gradient selection is represented analytically by coherence-order selection rather than explicit gradient pulses. The implementation is based on the Bruker `hsqcedetgp` pulse program and standard HSQC sequence.

## Physical / mathematical content

The sequence starts from longitudinal magnetisation on the F2 spin and uses two INEPT periods of `abs(1/(2*J))/2`, separated by simultaneous F1/F2 refocusing pulses. After the proton trim and transfer pulses, it evolves F1 in two halves with configured refocusing pulses and selects opposite F1 coherence orders for the echo and antiecho pathways. Each pathway then undergoes two multiplicity-editing intervals of `edit_time`, separated by another simultaneous refocusing pulse, followed by back-transfer, F2 coherence selection, optional F2 decoupling, and detection.

## Numerical / algorithmic content

The code forms `L = H + 1i*R + 1i*K`, uses dwell times `1/sweep(1)` and `1/sweep(2)`, and delegates propagation, pulses, and coherence selection to Spinach's `evolution`, `step`, and `coherence` routines. The function requires the `sphten-liouv` formalism.

## Syntax

```matlab
fid=hsqcedetgp(spin_system,parameters,H,R,K)
```

## Parameters / inputs

- `parameters.sweep`: two positive real sweep widths, `[F1 F2]`, Hz.
- `parameters.npoints`: two positive integer point counts, `[F1 F2]`.
- `parameters.spins`: two different isotope strings, `{F1 F2}`, present in the spin system (e.g. `{'13C','1H'}`).
- `parameters.decouple_f2`: cell array of isotopes to decouple in F2 (e.g. `{'15N','13C'}`).
- `parameters.decouple_f1`: cell array of isotopes receiving midpoint 180-degree refocusing pulses in F1 (e.g. `{'1H','15N'}`); it must not include the active F1 isotope.
- `parameters.J`: non-zero real scalar working scalar coupling, Hz.
- `parameters.trim_angle`: finite real proton trim-pulse angle, radians.
- `parameters.edit_time`: positive real multiplicity-editing delay, seconds; applied in two intervals.
- `H`: Hamiltonian matrix.
- `R`: relaxation superoperator.
- `K`: kinetics superoperator. These matrices must have matching dimensions.

## Outputs

- `fid.pos` and `fid.neg`: echo and antiecho signal components.

For natural-abundance simulations, the source recommends isotope dilution; see [`dilute.m`](https://spindynamics.org/wiki/index.php?title=dilute.m).

## References and links

- [HSQC sequence reference](https://doi.org/10.1016/0009-2614(80)80041-8)
- [HSQC reference](https://doi.org/10.1002/cmr.a.10095)
- [Multiplicity-edited HSQC reference](https://doi.org/10.1002/mrc.1260310315)
- [Spinach documentation for `hsqcedetgp.m`](https://spindynamics.org/wiki/index.php?title=hsqcedetgp.m)
