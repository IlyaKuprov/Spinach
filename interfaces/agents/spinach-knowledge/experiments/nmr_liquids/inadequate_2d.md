# experiments/nmr_liquids/inadequate_2d.m

- Signature: `fid=inadequate_2d(spin_system,parameters,H,R,K)`

## Purpose

Two-dimensional INADEQUATE. The indirect dimension encodes double-quantum frequencies and the direct dimension detects the converted signal. References: [the first cited paper](https://doi.org/10.1021/ja00398a044) and [the second cited paper](https://doi.org/10.1016/0022-2364(81)90060-3).

## Sequence and signal

The implementation starts from longitudinal magnetisation, uses J-coupling delays and pulses to create double-quantum coherence, selects the +2 and -2 coherence orders, and evolves them in F1 before F2 acquisition. Cosine and sine States pathways are retained. The delay is `abs(1/(4*parameters.J))`; decoupling is optional and defaults to an empty set.

## Syntax

```matlab
fid=inadequate_2d(spin_system,parameters,H,R,K)
```

## Parameters / inputs

- `parameters.sweep`: [F1 F2] sweep widths, Hz.
- `parameters.npoints`: [F1 F2] numbers of points.
- `parameters.spins`: active nucleus, e.g. `{'13C'}`.
- `parameters.decouple`: optional nuclei to decouple, e.g. `{'1H'}`; defaults to `{}`.
- `parameters.J`: working scalar coupling, Hz.
- `H`, `R`, `K`: Hamiltonian matrix, relaxation superoperator, and kinetics superoperator supplied by the context function.

## Outputs

- `fid.cos`, `fid.sin`: cosine and sine components of the States signal, with a double-quantum frequency axis in F1.

[Spin Dynamics Wiki page](https://spindynamics.org/wiki/index.php?title=inadequate_2d.m).
