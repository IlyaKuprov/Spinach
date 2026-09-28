# experiments/acquire.m

- Signature: `fid=acquire(spin_system,parameters,H,R,K)`

## Purpose

Acquire a free-induction decay by evolving the initial state under the supplied Hamiltonian, relaxation, and kinetics matrices, then observing it with the coil state. An optional dead time evolves the state before acquisition; optional homodecoupling modifies the Liouvillian during the sequence.

## Physical / mathematical content

The function forms `L = H + 1i*R + 1i*K`. It optionally applies analytical decoupling, adds `2*pi*homodec_pwr*homodec_oper` when homodecoupling is configured, and optionally advances the initial state for `dead_time`. The resulting signal is the observable evolution of `rho0` with `coil` as detection state.

## Numerical / algorithmic content

The acquisition interval is `1/sweep` seconds, and the function requests `npoints - 1` evolution steps after any optional dead-time evolution. The function delegates propagation to Spinach's `step` and `evolution` routines.

## Syntax

```matlab
fid=acquire(spin_system,parameters,H,R,K)
```

## Parameters / inputs

- `parameters.sweep`: positive real sweep width, Hz.
- `parameters.npoints`: positive integer number of FID points.
- `parameters.rho0`: initial state.
- `parameters.coil`: detection state.
- `parameters.decouple`: cell array of isotope strings to decouple; may be empty (for example, `{'15N','13C'}`). Non-empty analytical decoupling is limited to the `sphten-liouv` formalism.
- `parameters.homodec_oper`: optional square numeric operator added at the detection stage; if supplied, `parameters.homodec_pwr` is also required.
- `parameters.homodec_pwr`: optional real scalar power coefficient, Hz; must be paired with `homodec_oper`.
- `parameters.dead_time`: optional non-negative real delay, seconds, evolved before acquisition.
- `H`: Hamiltonian matrix.
- `R`: relaxation superoperator.
- `K`: kinetics superoperator. The three matrices must have matching dimensions.

## Outputs

- `fid`: FID observed in the state specified by `parameters.coil`.

## References and links

- [Spinach documentation for `acquire.m`](https://spindynamics.org/wiki/index.php?title=acquire.m)
