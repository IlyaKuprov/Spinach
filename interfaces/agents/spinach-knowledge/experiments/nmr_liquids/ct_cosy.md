# experiments/nmr_liquids/ct_cosy.m

- Signature: `fid=ct_cosy(spin_system,parameters,H,R,K)`

## Purpose

Constant-time COSY pulse sequence with analytical coherence selection, as described in [DOI 10.1006/jmra.1994.1095](https://doi.org/10.1006/jmra.1994.1095), [DOI 10.1080/00387010009350054](https://doi.org/10.1080/00387010009350054), and [DOI 10.1002/chem.201406283](https://doi.org/10.1002/chem.201406283). The source reverses the F1 trace at the end to preserve the conventional indirect-dimension sign for its reversed-delay implementation.

## Physical / mathematical content

- The sequence applies a 90-degree x pulse, selects +1 coherence, then samples a constant-time echo: for each t1 it evolves for `CT/2 - t1/2`, applies an x-axis pi pulse, and evolves for `CT/2 + t1/2`, with `CT` equal to the maximum t1.
- A final x pulse with angle `parameters.angle` precedes F2 acquisition. The initial state defaults to longitudinal magnetisation on `parameters.spins{1}`; the detection operator defaults to its `L+` state.

## Numerical / algorithmic content

- The F1 grid is `(0:npoints(1)-1)/sweep(1)` and the F2 dwell is `1/sweep(2)`. The source forms `L = H + 1i*R + 1i*K`, uses a `parfor` loop for the F1 values, and flips the acquired FID columns with `fliplr`. It requires the `sphten-liouv` formalism.

## Parameters / inputs

- `parameters.sweep`: two sweep widths `[F1 F2]` in Hz.
- `parameters.npoints`: two point counts `[F1 F2]`.
- `parameters.spins`: spin label used by the sequence, e.g. `'1H'` or `'13C'`.
- `parameters.angle`: final pulse angle in radians; defaults to `pi/2`.
- Optional `parameters.rho0` and `parameters.coil` override the default initial and detection states.
- `H`, `R`, and `K`: Hamiltonian, relaxation, and kinetics matrices received from the context function.

## Outputs

- `fid`: two-dimensional free induction decay, with the F1 direction flipped as described above.

[Spinach Wiki: ct_cosy.m](https://spindynamics.org/wiki/index.php?title=ct_cosy.m)
