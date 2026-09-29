# examples/fundamentals/mas_fplanck_powder.m

- Signature: `mas_fplanck_powder()`
- Source: [`examples/fundamentals/mas_fplanck_powder.m`](https://github.com/IlyaKuprov/Spinach/blob/main/examples/fundamentals/mas_fplanck_powder.m)

## Purpose

Compares a powder-averaged `13C` CSA signal under transverse RF between midpoint-sliced MAS (`rotor_stack`) and Liouville-space Fokker–Planck MAS (`singlerot`). Both routes use the same weighted two-angle Lebedev crystallite grid. The sliced route averages 7 and 13 rotor start phases explicitly; the Fokker–Planck route uses its uniform powder rotor-phase initial state.

## Convergence checks

The RF-assisted `L+` return signal is normalised by its initial overlap. Fokker–Planck rotor ranks 6 and 8, midpoint slice counts 33 and 65, and explicit phase grids 7 and 13 are refined independently; the `leb_2ang_rank_5` powder grid is compared with `leb_2ang_rank_11`. The sliced route visits rotor phases in the Fokker–Planck direction. Its entire `rotor_stack` spans one rotor period, so evolution lasts exactly one period. The test throws unless these refinement increments and the cross-route signal gap are below its stated targets, and prints `MAS_FPLANCK_POWDER_SUCCESS` on success.

This tests a finite-grid powder average of CSA, not a full quadrupolar powder satellite pattern or optimal-control gradients. See [`mas_fplanck_slices.m`](mas_fplanck_slices.md) for phase-sensitive single-crystal CSA and third-order quadrupolar central-transition checks.
