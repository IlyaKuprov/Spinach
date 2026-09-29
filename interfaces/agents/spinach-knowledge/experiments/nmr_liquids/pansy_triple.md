# experiments/nmr_liquids/pansy_triple.m

- Signature: `fid=pansy_triple(spin_system,parameters,H,R,K)`

## Purpose and pathway

Triple-channel PANSY, in magnitude-mode analytical-pathway form. The source starts from `Lz` on working spin 1, forms `L = H + 1i*R + 1i*K`, and creates first-pulse branches with `+pi/2` and `-pi/2` rotations about spin 1's `Lx`. It selects `+1` coherence on spin 1 and evolves both branches during F1. The second `pi/2` pulse uses the sum of the `Ly` operators for all three spins; branch subtraction is the stated axial-peak elimination. The source does not apply a later coherence-order filter. It notes that the working nuclei cannot be decoupled; gradient echo/anti-echo variants are separate sequences.

## Inputs and acquisition

- `parameters.spins`: three working nuclei in a cell array; source example: `{'1H','13C','15N'}`.
- `parameters.sweep`: three positive sweep widths in Hz; `parameters.npoints`: three positive integer point counts, ordered F1 then detection channels.
- `H`, `R`, and `K`: numeric, same-sized matrices supplied by the context function. The source requires the `sphten-liouv` formalism.

Each FID field is a two-dimensional array with `npoints(1)` F1-state columns: `fid.aa` detects `L+` on spin 1 with `npoints(1)` samples (F1/F1); `fid.ab` detects on spin 2 with `npoints(2)` samples (F1/F2); and `fid.ac` detects on spin 3 with `npoints(3)` samples (F1/F3). Their shapes are respectively `npoints(1) × npoints(1)`, `npoints(2) × npoints(1)`, and `npoints(3) × npoints(1)` (direct-dimension samples in rows).

## References

- [PANSY paper](https://doi.org/10.1021/ja0634876)
- [PANSY chapter](https://doi.org/10.1007/128_2011_226)
- [Canonical MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/nmr_liquids/pansy_triple.m)
- [Spinach Wiki: pansy_triple.m](https://spindynamics.org/wiki/index.php?title=pansy_triple.m)
