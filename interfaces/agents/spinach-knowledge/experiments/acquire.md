# experiments/acquire.m

- Source: [experiments/acquire.m](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/acquire.m)
- Signature: `fid=acquire(spin_system,parameters,H,R,K)`

## Purpose

Acquire a free-induction decay from a caller-supplied Spinach model. This routine consumes the initial state, detection state, Hamiltonian, relaxation superoperator, and kinetics superoperator; it does not prescribe a particular radical-pair/CIDNP mechanism or construct those physical inputs.

## Propagation and signal

The routine composes `L = H + 1i*R + 1i*K`, applies configured analytical decoupling to `L` and `parameters.rho0`, and optionally adds the homodecoupling term `2*pi*homodec_pwr*homodec_oper` (after projecting the operator into the Fokker–Planck space). If `parameters.dead_time` is present, it first propagates the initial state for that many seconds under the composed Liouvillian. It then observes the evolving state with `parameters.coil`, using an acquisition interval of `1/parameters.sweep` seconds and `parameters.npoints - 1` evolution steps.

## Parameters / inputs

- `parameters.rho0`: caller-provided initial state.
- `parameters.coil`: caller-provided detection state.
- `parameters.sweep`: sweep width in Hz.
- `parameters.npoints`: number of FID points.
- `parameters.decouple`: nuclei to decouple, for example `{'15N','13C'}`.
- Optional `parameters.homodec_oper` and `parameters.homodec_pwr`: operator and power coefficient; the source documents the power in Hz.
- Optional `parameters.dead_time`: pre-acquisition evolution time in seconds.
- `H`, `R`, and `K`: Hamiltonian, relaxation, and kinetics matrices supplied by the context function; these matrices must have matching dimensions.

The source defines a signal-acquisition operation, not a parameter sweep, a singlet-yield observable, or the chemistry encoded in `H`, `R`, and `K`.

## Output

- `fid`: FID observed in the state specified by `parameters.coil`.

## References and links

- [Spinach documentation for acquire.m](https://spindynamics.org/wiki/index.php?title=acquire.m)
