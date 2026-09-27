# experiments/nmr_liquids/inept.m

- Signature: `fid=inept(spin_system,parameters,H,R,K)`

## Purpose

Non-refocused INEPT. This variant returns the directly acquired coupled antiphase spectrum; it is not the refocused, broadband-decoupled INEPT variant. The source cites [this paper](https://doi.org/10.1021/ja00497a058).

## Sequence and signal

The sequence starts from isotropic thermal equilibrium, applies pulses to the two working spin channels with J-coupling evolution intervals of `abs(1/(4*parameters.J))`, and uses phase-cycled pulses on the second configured spin channel before direct acquisition. The detected result is the coupled antiphase FID, rather than a refocused broadband-decoupled spectrum.

## Syntax

```matlab
fid=inept(spin_system,parameters,H,R,K)
```

## Parameters / inputs

- `parameters.sweep`: F1 sweep width, Hz.
- `parameters.npoints`: number of points.
- `parameters.spins`: {F1 F2} working nuclei, e.g. `{'15N','1H'}`.
- `parameters.J`: working scalar coupling, Hz.
- `H`, `R`, `K`: Hamiltonian matrix, relaxation superoperator, and kinetics superoperator supplied by the context function.

## Output

- `fid`: directly acquired free induction decay.

[Spin Dynamics Wiki page](https://spindynamics.org/wiki/index.php?title=inept.m).
