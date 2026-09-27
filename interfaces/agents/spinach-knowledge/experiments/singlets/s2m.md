# experiments/singlets/s2m.m

- Signature: `rho=s2m(spin_system,L,Hx,Hy,rho,J,delta_v)`

## Purpose

Implements the S2M sequence of Pileio and Levitt.

## Numerical / algorithmic content

The routine sets `t=1/(4*sqrt(J^2+delta_v^2))` and rounds `floor(pi*abs(J)/(2*abs(delta_v)))` up to an even repetition count `m1`. It first runs `m1/2` blocks of evolution under `L` for `t`, a `pi` pulse about `Hx`, and another `L` evolution for `t`. It then evolves under `L` for `t`, applies a `sign(J)*pi/2` pulse about `Hx`, runs `m1` such evolution/pulse/evolution blocks, and finishes with a `pi/2` pulse about `Hy`.

## Parameters / inputs

- `L` — background Liouvillian
- `Hx` — X spin operator
- `Hy` — Y spin operator
- `rho` — initial state vector
- `J` — J-coupling (Hz); the phase of the 90-degree pulse next to the lone tau delay follows its sign
- `delta_v` — Zeeman frequency difference (Hz)

## Outputs

- `rho` — final state vector

## Reference

- [Spinach documentation](https://spindynamics.org/wiki/index.php?title=s2m.m)
