# experiments/nmr_liquids/cosy.m

- Signature: `fid=cosy(spin_system,parameters,H,R,K)`

## Purpose

Phase-sensitive COSY sequence; the source cites [DOI 10.1063/1.432450](https://doi.org/10.1063/1.432450) and [DOI 10.1016/0022-2364(82)90279-7](https://doi.org/10.1016/0022-2364(82)90279-7).

## Physical / mathematical content

- The sequence prepares longitudinal magnetisation on `parameters.spins{1}`, applies a 90-degree x pulse, evolves along F1, and selects the +1 coherence pathway. It then applies the second x pulse with angle `parameters.angle` and acquires F2 with the same spin as the detected observable.

## Numerical / algorithmic content

- Both dimensions use the reciprocal of the scalar sweep width. The F1 trajectory has `parameters.npoints(1)` points; F2 acquisition has `parameters.npoints(2)` points. The source combines the inputs as `L = H + 1i*R + 1i*K` and requires the `sphten-liouv` formalism.

## Syntax

```matlab
fid=cosy(spin_system,parameters,H,R,K)
```

## Parameters / inputs

- `parameters.sweep`: scalar sweep width in Hz.
- `parameters.npoints`: point counts for F1 and F2.
- `parameters.spins`: spin label used by the sequence, e.g. `'1H'` or `'13C'`.
- `parameters.angle`: second-pulse angle in radians; the source notes COSY45 and COSY60 variants as examples.
- `H`, `R`, and `K`: Hamiltonian, relaxation, and kinetics matrices received from the context function.

## Outputs

- `fid`: two-dimensional free induction decay.
- The implementation analytically retains the +1 t1 coherence order, corresponding to one phase-sensitive pathway. If the second pulse is not 90 degrees, magnitude-mode plotting is advised.

[Spinach Wiki: cosy.m](https://spindynamics.org/wiki/index.php?title=cosy.m)
