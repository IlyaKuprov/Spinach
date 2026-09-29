# experiments/nmr_liquids/inadequate_2d.m

- Signature: `fid=inadequate_2d(spin_system,parameters,H,R,K)`

## Purpose

Two-dimensional INADEQUATE. The source describes the F1 axis as a double-quantum frequency coordinate and cites [the first paper](https://doi.org/10.1021/ja00398a044) and [the second paper](https://doi.org/10.1016/0022-2364(81)90060-3). It recommends `dilute.m` for generating carbon-pair isotopomers.

## Sequence and signal

The source initialises longitudinal magnetisation on the configured nucleus and uses two delays `tau=abs(1/(4*parameters.J))` around an inversion pulse. A final 90-degree pulse creates the double-quantum signal, which is filtered to coherence orders +2 and -2. The F1 evolution is sampled at dwell `1/parameters.sweep(1)`; a following 90-degree pulse converts the stored states for detection, and F2 is acquired at dwell `1/parameters.sweep(2)`. Sweep widths and J are in Hz, so `tau` is in seconds. The returned `fid.cos` and `fid.sin` are the two States quadrature components; each is a two-dimensional array with `npoints(2)` F2 acquisition samples in rows and `npoints(1)` F1 increments in columns; the source passes the F1 state stack directly to `evolution` without transposing the result.

## Inputs

- `parameters.sweep`: two positive sweep widths, `[F1 F2]`, in Hz.
- `parameters.npoints`: two positive integer point counts, `[F1 F2]`.
- `parameters.spins`: active nucleus in a cell array; source example: `{'13C'}`.
- `parameters.decouple`: optional cell array of nuclei to decouple; source example: `{'1H'}`. If absent, the function sets it to an empty cell array.
- `parameters.J`: working scalar coupling in Hz.
- `H`, `R`, `K`: Hamiltonian matrix, relaxation superoperator, and kinetics superoperator supplied by the context function. The source requires the `sphten-liouv` formalism and same-sized matrices.

[Canonical MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/nmr_liquids/inadequate_2d.m) · [Spin Dynamics Wiki page](https://spindynamics.org/wiki/index.php?title=inadequate_2d.m).
