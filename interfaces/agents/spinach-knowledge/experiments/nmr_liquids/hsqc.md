# experiments/nmr_liquids/hsqc.m

- Signature: `fid=hsqc(spin_system,parameters,H,R,K)`

## Purpose

Phase-sensitive heteronuclear single-quantum coherence (HSQC). The sequence uses scalar-coupling evolution to transfer coherence between the two configured spin channels, records the States quadrature signal, and acquires F2. References: [HSQC paper](https://doi.org/10.1016/0009-2614(80)80041-8) and [review](https://doi.org/10.1002/cmr.a.10095).

## Sequence and signal

The coupling delay is set from `parameters.J` as `abs(1/(2*J))`. The sequence applies the F2 excitation, J-coupling evolution and inversion pulses, then constructs the indirect F1 evolution. `parameters.decouple_f1` specifies nuclei receiving the midpoint 180-degree refocusing pulses in F1; `parameters.decouple_f2` specifies nuclei decoupled during F2. The two States quadrature pathways are returned separately.

If omitted, `parameters.rho0` is initialized to an F2 longitudinal state and `parameters.coil` to an F2 raising-operator detection state. Natural-abundance simulations should use isotope dilution; see `dilute.m`.

## Syntax

```matlab
fid=hsqc(spin_system,parameters,H,R,K)
```

## Parameters / inputs

- `parameters.sweep`: [F1 F2] sweep widths, Hz.
- `parameters.npoints`: [F1 F2] numbers of points.
- `parameters.spins`: {F1 F2} nuclei, e.g. `{'13C','1H'}`.
- `parameters.decouple_f2`: nuclei to decouple in F2, e.g. `{'15N','13C'}`.
- `parameters.decouple_f1`: nuclei receiving midpoint 180-degree refocusing pulses in F1, e.g. `{'1H','13C'}`.
- `parameters.J`: working scalar coupling, Hz.
- Optional `parameters.rho0` and `parameters.coil`: initial state and detection state.
- `H`, `R`, `K`: Hamiltonian matrix, relaxation superoperator, and kinetics superoperator supplied by the context function.

## Outputs

- `fid.pos`, `fid.neg`: the two components of the States quadrature signal.

[Spin Dynamics Wiki page](https://spindynamics.org/wiki/index.php?title=hsqc.m).
