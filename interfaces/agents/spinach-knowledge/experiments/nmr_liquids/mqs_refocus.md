# experiments/nmr_liquids/mqs_refocus.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/nmr_liquids/mqs_refocus.m) · [Spin Dynamics Wiki](https://spindynamics.org/wiki/index.php?title=mqs_refocus.m)

- Signature: `fid=mqs_refocus(spin_system,parameters,H,R,K)`

## Purpose and sequence

A two-dimensional multiple-quantum correlation sequence with refocusing. The source cites [the original article](https://doi.org/10.1063/1.432450) and [the earlier sequence paper](https://doi.org/10.1016/0022-2364(80)90096-7). It is homonuclear in practice and uses exact analytical coherence-order projection rather than an explicit phase cycle.

Starting from `rho0`, it applies a 90° x pulse, evolves for `delay_1`, applies a 180° x pulse, and repeats `delay_1`. The next 90° pulse is x for even `mqorder(1)` and y for odd; the source selects this coherence order on `spins{1}`. After F1 trajectory evolution and the user-set x pulse `angle`, it selects `mqorder(2)` on the same spin, evolves for `delay_2`, applies a 180° x refocusing pulse, repeats `delay_2`, and acquires F2.

## Inputs

- `parameters.sweep`: two positive sweep widths in Hz, ordered [F1 F2]; `parameters.npoints`: two positive integer point counts in that order.
- `parameters.spins`: two labels; the header example is `{'1H','1H'}`.
- `parameters.angle`: flip angle in radians; `parameters.mqorder`: the two coherence orders selected at the two explicit projection steps.
- `parameters.delay_1` and `parameters.delay_2`: sequence delays in seconds.
- `parameters.rho0`: initial state; `parameters.coil`: detection state.
- `H`, `R`, and `K`: Hamiltonian, relaxation, and kinetics matrices. The source validates equal matrix dimensions and requires `sphten-liouv` formalism.

## Output

`fid` is a two-dimensional free induction decay for amplitude-mode processing. The F1 trajectory is stacked as initial states and F2 is the observable evolution, giving array shape `[npoints(2), npoints(1)]` (F2 rows, F1 columns).
