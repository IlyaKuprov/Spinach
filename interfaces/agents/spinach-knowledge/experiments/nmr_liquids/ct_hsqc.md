# experiments/nmr_liquids/ct_hsqc.m

- Signature: `fid=ct_hsqc(spin_system,parameters,H,R,K)`

## Purpose

Constant-time phase-sensitive HSQC sequence, citing [DOI 10.1016/0022-2364(92)90144-V](https://doi.org/10.1016/0022-2364(92)90144-V) and [DOI 10.1007/BF00227470](https://doi.org/10.1007/BF00227470).

## Physical / mathematical content

- The source sets the J-evolution interval to `abs(1/(4*parameters.J))`, prepares the initial state on spin 2, and applies the transfer and refocusing pulses on the two specified spins. It samples a constant-time t1 grid and separates the two States quadrature pathways as `fid.pos` and `fid.neg`.
- The requested F2 decoupling is applied before detection; the detection state defaults to `L+` on `parameters.spins{2}`. The implementation forms `L = H + 1i*R + 1i*K`; the initial state defaults to longitudinal magnetisation on `parameters.spins{2}`, and the detection state to its `L+` operator.

## Numerical / algorithmic content

- The two sweep widths set the F1 grid and F2 dwell (`1/sweep(2)`). The source requires the `sphten-liouv` formalism and equal-sized matrix inputs `H`, `R`, and `K`; it does not perform orientation or geometry averaging.

## Syntax

```matlab
fid=ct_hsqc(spin_system,parameters,H,R,K)
```

## Parameters / inputs

- `parameters.sweep`: two sweep widths `[F1 F2]` in Hz.
- `parameters.npoints`: two point counts `[F1 F2]`.
- `parameters.spins`: two spin labels `{F1 F2}`, e.g. `'13C'` and `'1H'`.
- `parameters.decouple_f2`: optional nuclei to decouple during F2; defaults to empty (for example, `{'15N','13C'}`).
- `parameters.J`: working scalar coupling in Hz.
- Optional `parameters.rho0` and `parameters.coil` set the initial and detection states.
- `H`, `R`, and `K`: Hamiltonian, relaxation, and kinetics matrices received from the context function.

## Outputs

- `fid.pos` and `fid.neg`: the two components of the States quadrature signal. For natural-abundance simulations, use isotope dilution; see `dilute.m`.

[Spinach Wiki: ct_hsqc.m](https://spindynamics.org/wiki/index.php?title=ct_hsqc.m)
