# experiments/bruker/hmqcetgpsi.m

- Signature: `fid=hmqcetgpsi(spin_system,parameters,H,R,K)`

## Purpose

Simulate the sensitivity-improved echo/antiecho, gradient-selected HMQC sequence. Gradient selection is represented analytically by coherence-order selection rather than explicit gradient pulses. The implementation follows the Bruker `hmqcetgpsi` pulse program and standard HMQC sequence.

## Physical / mathematical content

The sequence begins with longitudinal magnetisation on the F2 spin and a J-coupling transfer period `abs(1/(2*J))`. It acquires the indirect dimension as two evolution halves with configured F1 refocusing pulses, then selects opposite F1 coherence orders for the two pathways. The sensitivity-improvement section adds three further J-evolution periods with proton/carbon pulse pairs and a final proton echo pulse before F2 coherence selection and detection.

## Numerical / algorithmic content

The code forms `L = H + 1i*R + 1i*K`, uses dwell times `1/sweep(1)` and `1/sweep(2)`, and delegates propagation and pulses to Spinach's `evolution` and `step` routines. The function requires the `sphten-liouv` formalism.

## Syntax

```matlab
fid=hmqcetgpsi(spin_system,parameters,H,R,K)
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
- [Sensitivity-improvement sequence reference](https://doi.org/10.1016/0022-2364(91)90036-S)
- [Spinach documentation for `hmqcetgpsi.m`](https://spindynamics.org/wiki/index.php?title=hmqcetgpsi.m)
