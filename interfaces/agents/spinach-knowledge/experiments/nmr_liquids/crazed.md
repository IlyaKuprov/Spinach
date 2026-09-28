# experiments/nmr_liquids/crazed.m

- Signature: `fid=crazed(spin_system,parameters,H,R,K)`

## Purpose

CRAZED pulse sequence, implemented as an ideal analytical coherence-pathway version of the sequence described in [the cited paper](https://doi.org/10.1126/science.8266096). The gradient selection is represented by explicit coherence projections.

## Physical / mathematical content

- Starting from `parameters.rho0`, the sequence applies a 90-degree y pulse to `parameters.spins{1}`, evolves through the F1 trajectory, and selects coherence order +2 (the double-quantum branch).
- It applies the second y pulse with angle `parameters.angle`, selects coherence order +1 (the observable single-quantum branch), then evolves and detects on the same spin.

## Numerical / algorithmic content

- The scalar time step is `1/parameters.sweep`. The F1 trajectory contains `parameters.npoints(1)` points and the F2 acquisition contains `parameters.npoints(2)` points. The source forms `L = H + 1i*R + 1i*K` and requires the `sphten-liouv` formalism.

## Syntax

```matlab
fid=crazed(spin_system,parameters,H,R,K)
```

## Parameters / inputs

- `parameters.sweep`: sweep width in Hz.
- `parameters.npoints`: point counts for both dimensions.
- `parameters.spins`: spin label used by the sequence, e.g. `'1H'` or `'13C'`.
- `parameters.angle`: second pulse angle.
- `parameters.rho0`: initial condition.
- `H`, `R`, and `K`: Hamiltonian, relaxation, and kinetics matrices received from the context function.

## Outputs

- `fid`: two-dimensional free induction decay. The source models gradient selection with explicit projection onto the +2 double-quantum branch followed by the +1 observable single-quantum branch.

[Spinach Wiki: crazed.m](https://spindynamics.org/wiki/index.php?title=crazed.m)
