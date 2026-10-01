# examples/fundamentals/mas_fplanck_slices.m

- Signature: `mas_fplanck_slices()`
- Source: [`examples/fundamentals/mas_fplanck_slices.m`](https://github.com/IlyaKuprov/Spinach/blob/main/examples/fundamentals/mas_fplanck_slices.m)

## Purpose

Tests agreement between explicitly midpoint-sliced MAS propagation (`rotor_stack`) and the Liouville-space Fokker–Planck rotor generator (`singlerot`). Both routes start from the same single crystal at phase zero, use the same magic-angle rotor and spin Hamiltonian, and independently refine their discretisations.

## Cases and method

The `13C` case uses an L+ to Lz transfer under anisotropic shielding and continuous transverse RF; this phase-sensitive observable distinguishes opposite rotor traversal directions. The `27Al` case uses the 3.2 MHz, asymmetry-0.16 quadrupolar tensor and shielding of the Šmelko example, with the central-transition coherence selectively prepared and the rotating-frame transformation taken through third order. Its third-order term is checked against the second-order generator to ensure it is nonzero.

A single normalised complex signal is evaluated after one rotor period. The FP computation applies a phase-point delta and sums the detection operator over rotor phases. For a positive rotor rate, the FP phase delta travels toward decreasing phase. The sliced computation therefore starts at the negative half-step and visits rotor-stack phases in reverse order; each route includes the same carrier transformation. It reports rank/slice signals, last refinement increments, and their difference, and throws if either refinement or the cross-route difference exceeds its normalised-signal target (0.001 for CSA, 0.0001 for the quadrupolar central transition).

## Output and limits

The script prints `MAS_FPLANCK_SLICES_SUCCESS` only after both assertions and the third-order-term check pass; it produces no figure. These tests establish equivalence for the specified single-crystal observables and settings, not for a powder average, the entire strongly quadrupolar satellite manifold, or a complete optimal-control gradient. The separate [`mas_fplanck_powder.m`](mas_fplanck_powder.md) checks powder-averaged CSA on a weighted Lebedev grid.
