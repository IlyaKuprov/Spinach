# experiments/singlets/m2s.m

- Signature: `rho=m2s(spin_system,L,Hx,Hy,rho,J,delta_v)`

## Purpose

Implements the M2S sequence of Pileio and Levitt.

## Numerical / algorithmic content

The routine sets `t=1/(4*sqrt(J^2+delta_v^2))` and chooses an even repetition count from `floor(pi*abs(J)/(2*abs(delta_v)))`. It applies a `pi/2` pulse about `Hy`, alternates evolution under `L` for `t` with `pi` pulses about `Hx` for the full repetition count, then applies a `sign(J)*pi/2` pulse about `Hx` and an additional `L` evolution. A final loop runs half as many evolution/`Hx`-pulse/evolution blocks.

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

- [Spinach documentation](https://spindynamics.org/wiki/index.php?title=m2s.m)
