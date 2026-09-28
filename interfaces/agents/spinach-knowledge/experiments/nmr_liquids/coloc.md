# experiments/nmr_liquids/coloc.m

- Signature: `fid=coloc(spin_system,parameters,H,R,K)`

## Purpose

COLOC NMR pulse sequence implementing Fig. 1b of [the cited paper](https://doi.org/10.1016/0022-2364(84)90136-7), with the dashed pulses during the delta(2) period omitted. The delay delta(1) defaults to half of the maximum F1 evolution time implied by the sweep width; delta(2) is specified by the caller and is typically `40e-3` seconds.

## Physical / mathematical content

- The initial state is longitudinal magnetisation on `parameters.spins{1}`; the sequence applies a 90-degree x pulse, then samples an embedded echo over the F1 grid. At each point it evolves for half the current t1, applies pi pulses about x to spins 1 and 2, then evolves for `delta1 - t1/2`.
- The resulting stack is selected for coherence order -1 on spin 1, rotated by 90 degrees about x on spin 1 and about y on spin 2, evolved for delta(2), and selected for coherence order +1 on spin 2. Spin 1 is decoupled and the signal is detected on spin 2.
- The effective Liouvillian is `L = H + 1i*R + 1i*K`, matching the function's Hamiltonian, relaxation, and kinetics inputs.

## Numerical / algorithmic content

- The indirect evolution times are `t1 = (0:npoints(1)-1)/sweep(1)`; the echo is split around the pi-pulse pair as described above. Direct-dimension acquisition uses the reciprocal second sweep width and `npoints(2)` points.
- The implementation requires the `sphten-liouv` formalism. `H`, `R`, and `K` must be numeric matrices of the same dimensions.

## Parameters / inputs

- `parameters.sweep`: two positive sweep widths `[F1 F2]`, in Hz.
- `parameters.npoints`: two point counts `[F1 F2]`.
- `parameters.spins`: two spin labels `{F1 F2}` (for example, `'13C'` and `'1H'`).
- `parameters.delta2`: COLOC delta(2) delay, typically `40e-3` seconds.
- `parameters.delta1`: optional COLOC delta(1) delay in seconds; when supplied it must be at least half the maximum t1. If omitted, the implementation sets it to `(npoints(1)-1)/(2*sweep(1))`.
- `H`, `R`, and `K`: same-sized Hamiltonian, relaxation, and kinetics matrices, respectively, supplied by the context function.

## Outputs

- `fid`: free induction decay for magnitude-mode processing.
- Natural-abundance simulations should use isotope dilution; see `dilute.m`.

[Spinach Wiki: coloc.m](https://spindynamics.org/wiki/index.php?title=coloc.m)
