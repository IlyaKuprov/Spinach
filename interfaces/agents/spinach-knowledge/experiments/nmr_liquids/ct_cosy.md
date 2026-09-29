# experiments/nmr_liquids/ct_cosy.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/nmr_liquids/ct_cosy.m) · [Spinach Wiki: ct_cosy.m](https://spindynamics.org/wiki/index.php?title=ct_cosy.m)

- Signature: `fid=ct_cosy(spin_system,parameters,H,R,K)`

## Purpose

Constant-time COSY with analytical coherence selection, described in [DOI 10.1006/jmra.1994.1095](https://doi.org/10.1006/jmra.1994.1095), [DOI 10.1080/00387010009350054](https://doi.org/10.1080/00387010009350054), and [DOI 10.1002/chem.201406283](https://doi.org/10.1002/chem.201406283). The implementation reverses the F1 trace at the end to preserve the conventional indirect-dimension sign for its reversed-delay implementation.

## Inputs

- `parameters.sweep`: two positive sweep widths in Hz, ordered [F1, F2].
- `parameters.npoints`: two positive integer point counts, ordered [F1, F2].
- `parameters.spins`: a one-element cell array naming the observed isotope, such as `{'1H'}` or `{'13C'}`.
- `parameters.angle`: final pulse angle in radians; if omitted, the source sets it to `pi/2`.
- `parameters.rho0`: optional initial state; if omitted, the source uses the `Lz` state of the selected spin.
- `parameters.coil`: optional detection state; if omitted, the source uses the `L+` state of the selected spin.
- `H`, `R`, and `K`: Hamiltonian, relaxation, and kinetics matrices supplied by the context function; they form `H + 1i*R + 1i*K`. The implementation requires the `sphten-liouv` formalism.

## Sequence and coherence selection

The sequence applies a 90-degree x pulse and selects +1 coherence. For each F1 point `t1`, it evolves for `CT/2 - t1/2`, applies a pi pulse, then evolves for `CT/2 + t1/2`, where `CT` is the last point on the source's F1 time grid. It applies the configurable final pulse and performs F2 observable evolution using timestep `1/parameters.sweep(2)`. The output is flipped left-to-right, reversing its F1 direction as noted above.

## Output and scope

- `fid`: two-dimensional free induction decay with F1 and F2 evolution; its F1 direction is reversed by the source's final `fliplr` operation.
