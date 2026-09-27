# experiments/bruker/hmqcetgp.m

- Signature: `fid=hmqcetgp(spin_system,parameters,H,R,K)`

## Purpose

Simulate the echo/antiecho, gradient-selected HMQC sequence. Gradient selection is represented analytically by coherence-order selection rather than explicit gradient pulses. The implementation is based on the Bruker `hmqcetgp` pulse program and standard HMQC sequence.

## Physical / mathematical content

The sequence starts from longitudinal magnetisation on the F2 spin, uses a J-coupling transfer period of `abs(1/(2*J))`, and evolves the indirect dimension in two halves separated by configured F1 refocusing pulses. It selects opposite F1 coherence orders for the echo and antiecho pathways, applies the back-transfer steps, selects F2 single-quantum coherence, and detects both signals with the F2 transverse state.

## Numerical / algorithmic content

The code constructs `L = H + 1i*R + 1i*K`, uses dwell times `1/sweep(1)` and `1/sweep(2)`, and delegates propagation and pulse actions to Spinach's `evolution` and `step` routines. The function requires the `sphten-liouv` formalism.

## Syntax

```matlab
fid=hmqcetgp(spin_system,parameters,H,R,K)
```

## Parameters / inputs

- `parameters.sweep`: two positive real sweep widths, `[F1 F2]`, Hz.
- `parameters.npoints`: two positive integer point counts, `[F1 F2]`.
- `parameters.spins`: two different isotope strings, `{F1 F2}`, present in the spin system (e.g. `{'13C','1H'}`).
- `parameters.decouple_f2`: cell array of isotopes to decouple in F2 (e.g. `{'15N','13C'}`).
- `parameters.decouple_f1`: cell array of isotopes receiving midpoint 180-degree refocusing pulses in F1 (e.g. `{'1H'}`); it must not include the active F1 isotope.
- `parameters.J`: non-zero real scalar working scalar coupling, Hz.
- `H`: Hamiltonian matrix.
- `R`: relaxation superoperator.
- `K`: kinetics superoperator. These matrices must have matching dimensions.

## Outputs

- `fid.pos` and `fid.neg`: echo and antiecho signal components.

For natural-abundance simulations, the source recommends isotope dilution; see [`dilute.m`](https://spindynamics.org/wiki/index.php?title=dilute.m).

## References and links

- [Standard HMQC sequence](https://doi.org/10.1016/0022-2364(83)90241-X)
- [Spinach documentation for `hmqcetgp.m`](https://spindynamics.org/wiki/index.php?title=hmqcetgp.m)
