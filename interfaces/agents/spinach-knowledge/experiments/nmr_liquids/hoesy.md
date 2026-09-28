# experiments/nmr_liquids/hoesy.m

- Signature: `fid=hoesy(spin_system,parameters,H,R,K)`

## Purpose

Phase-sensitive heteronuclear NOESY with indirect evolution in F1, a mixing period, and direct acquisition in F2. The implementation follows the sequence described in [the original paper](https://doi.org/10.1021/ja00353a071) and [this later reference](https://doi.org/10.1039/C8CP00911B).

## Sequence and signal

The function builds `L = H + 1i*R + 1i*K`, applies the first F1 pulse, and evolves the two indirect-time halves with the configured F1 decoupling and refocusing. Phase-cycled F1 pulses generate cosine and sine pathways; homospoil is applied before the mixing evolution, then an F2 pulse and direct acquisition produce the two States components. The mixing period propagates relaxation and kinetics (`1i*R + 1i*K`).

This is an ideal heteronuclear NOESY model: gradient and diffusion attenuation, finite-pulse losses, and experimental normalisation are not included.

## Syntax

```matlab
fid=hoesy(spin_system,parameters,H,R,K)
```

## Parameters / inputs

- `parameters.sweep`: two sweep widths [F1 F2], Hz.
- `parameters.npoints`: numbers of FID points [F1 F2].
- `parameters.spins`: nuclei used on F1 and F2, e.g. `{'15N','13C'}`.
- `parameters.decouple_f1`: nuclei to decouple in F1, e.g. `{'1H','13C'}`.
- `parameters.tmix`: mixing time, seconds.
- `parameters.rho0`: initial state.
- `parameters.needs`: set to `{'rho_eq'}`; the sequence requires the thermal-equilibrium state.
- `H`, `R`, `K`: Hamiltonian matrix, relaxation superoperator, and kinetics superoperator supplied by the context function.

## Outputs

- `fid.cos`, `fid.sin`: cosine and sine components of the F1 hypercomplex FID for States processing.

[Spin Dynamics Wiki page](https://spindynamics.org/wiki/index.php?title=hoesy.m).
