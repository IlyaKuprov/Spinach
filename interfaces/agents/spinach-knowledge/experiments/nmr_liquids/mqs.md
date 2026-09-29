# experiments/nmr_liquids/mqs.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/nmr_liquids/mqs.m) · [Spin Dynamics Wiki](https://spindynamics.org/wiki/index.php?title=mqs.m)

- Signature: `fid=mqs(spin_system,parameters,H,R,K)`

## Purpose and sequence

A two-dimensional multiple-quantum experiment in the non-refocused MQ/MaxQ variant. The source cites [this paper](https://doi.org/10.1002/cphc.201800667) and [this communication](https://doi.org/10.1039/d1cc03079e). It should be called through `liquid.m`, which supplies `H`, `R`, and `K`.

The function applies an initial 90° x pulse, evolves for `delay`, applies a 180° x pulse, and evolves for a second `delay`. The following 90° pulse is about x for an even `mqorder` and y for an odd one; then `coherence` explicitly retains that order on `spins{1}`. It records the first indirect dimension as a trajectory, applies the user-specified final flip angle about x, and detects the second dimension. The source does not apply a post-mixing refocusing period; use `mqs_refocus.m` for that variant.

## Inputs

- `parameters.sweep`: two positive sweep widths in Hz, ordered [F1 F2].
- `parameters.npoints`: two positive integer point counts, ordered [F1 F2].
- `parameters.spins`: two spin labels; the source-header example is `{'1H','1H'}`.
- `parameters.angle`: final flip angle in radians.
- `parameters.delay`: J-coupling evolution delay in seconds (used for both delays around the 180° pulse).
- `parameters.mqorder`: coherence order selected on the first listed spin.
- `parameters.rho0`: initial state; `parameters.coil`: detection state.
- `H`, `R`, and `K`: Hamiltonian, relaxation, and kinetics matrices supplied by the liquid context. The source requires matching matrix dimensions and the `sphten-liouv` formalism.

## Output

`fid` is a two-dimensional magnitude-mode free induction decay. The source requests F1 then F2 in `npoints`; the final `observable` evolution returns the direct F2 samples by rows and the F1 trajectory samples by columns, so the array shape is `[npoints(2), npoints(1)]`. The source does not report a measured spectrum or simulation result.
