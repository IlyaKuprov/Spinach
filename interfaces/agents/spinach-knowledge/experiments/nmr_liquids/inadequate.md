# experiments/nmr_liquids/inadequate.m

- Signature: `fid=inadequate(spin_system,parameters,H,R,K)`

## Purpose

INADEQUATE selects double-quantum coherence from coupled carbon pairs and converts it back into observable single-quantum magnetisation. At natural-abundance 13C, this gives 13C pair subspectra. The implementation is described in [the cited paper](https://doi.org/10.1021/ja00534a056); use `dilute.m` to generate carbon-pair isotopomers.

## Sequence and signal

The sequence uses J-coupling evolution and pulses to create and select double-quantum coherence, then converts the selected coherence back to detectable single-quantum magnetisation for the FID. The supplied Hamiltonian, relaxation, and kinetics operators are combined for propagation.

## Syntax

```matlab
fid=inadequate(spin_system,parameters,H,R,K)
```

## Parameters / inputs

- `parameters.sweep`: sweep width, Hz.
- `parameters.npoints`: number of FID points.
- `parameters.spins`: active nuclei, e.g. `{'13C'}`.
- `parameters.decouple`: nuclei to decouple, e.g. `{'1H'}`.
- `parameters.J`: working J-coupling, Hz.
- `H`, `R`, `K`: Hamiltonian matrix, relaxation superoperator, and kinetics superoperator supplied by the context function.

## Output

- `fid`: free induction decay.

[Spin Dynamics Wiki page](https://spindynamics.org/wiki/index.php?title=inadequate.m).
