# experiments/singlets/s2m.m

- MATLAB source: [experiments/singlets/s2m.m](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/singlets/s2m.m)
- Signature: `rho=s2m(spin_system,L,Hx,Hy,rho,J,delta_v)`

## What the routine does

The source identifies this as the S2M sequence of Pileio and Levitt. It accepts the caller's initial state `rho`; it does not construct a magnetisation preparation or apply an explicit coherence-order filter. The pulse and free-evolution steps are the sequence implementation, not a measurement of transfer efficiency or a promise of a particular final singlet population.

## Sequence and timing

The routine sets `t=1/(4*sqrt(J^2+delta_v^2))`. Because `J` and `delta_v` are in Hz, `t` is in seconds. It sets `m1=floor(pi*abs(J)/(2*abs(delta_v)))` and increments `m1` by one when it is odd, giving an even pulse-repeat count.

In order, it performs:

1. `m1/2` repetitions of free evolution under `L` for `t`, an `Hx` rotation by `pi`, and another `L` evolution for `t`.
2. One `L` evolution for `t`, followed by an `Hx` rotation by `sign(J)*pi/2`.
3. `m1` repetitions of `L` for `t`, an `Hx` rotation by `pi`, and `L` for `t`.
4. A final `Hy` rotation by `pi/2`.

The source uses `J` and `delta_v` (both in Hz) to set the delay and repetition count; their sign/magnitude are not inserted as additional Hamiltonian terms by this function. The sign of `J` selects the sign of the middle 90-degree X rotation. `L` is the supplied background Liouvillian; `Hx` and `Hy` are supplied spin-operator matrices.

## Inputs and output

- `L`, `Hx`, and `Hy` must be numeric matrices of matching dimensions; the row dimension of `rho` must match the column dimension of `L`.
- `J` and `delta_v` must each be finite, real, nonzero numeric scalars.
- Output `rho` is the final propagated state, with the state-space dimension of the input.

## References

- [Canonical MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/singlets/s2m.m)
- [Spinach Wiki: s2m.m](https://spindynamics.org/wiki/index.php?title=s2m.m)
