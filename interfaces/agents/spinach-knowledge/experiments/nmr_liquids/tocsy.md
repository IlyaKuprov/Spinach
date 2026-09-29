# experiments/nmr_liquids/tocsy.m

- Signature: `fid=tocsy(spin_system,parameters,H,R,K)`

## Purpose and sequence

Amplitude-mode homonuclear TOCSY with the source's continuous-spin-lock model. Starting from `parameters.rho0`, the routine applies a `pi/2` pulse about `Lx` on the working spin, then records F1 evolution under `L = H + 1i*R + 1i*K`. During the mixing interval, it propagates two branches with `L + 2*pi*lamp*Lx` and `L + 2*pi*lamp*Ly`; these include the full Liouvillian as well as the spin-lock term. It detects with `L+` on the same spin during F2 and returns `fid.cos` and `fid.sin`, the States-quadrature components, each with shape `npoints(2) × npoints(1)` (F2 direct-time samples in rows, F1 state-stack samples in columns).

This is not an explicit MLEV, DIPSI, WALTZ, or clean-TOCSY pulse-train simulation.

## Inputs

- `parameters.spins`: one homonuclear species in a cell array; source examples include `{'1H'}` and `{'13C'}`.
- `parameters.sweep`: two positive sweep widths in Hz; `parameters.npoints`: two positive integer point counts, in F1/F2 order.
- `parameters.tmix`: non-negative scalar mixing time in seconds; `parameters.lamp`: spin-lock power in Hz, used in the source as `2*pi*lamp` angular frequency; `parameters.rho0`: numeric initial state with a Liouville-space dimension matching `H`.
- `H`, `R`, and `K`: numeric, same-sized matrices from the context function. Supported formalisms are `sphten-liouv` and `zeeman-liouv`.

## References

- [TOCSY paper](https://doi.org/10.1016/0022-2364(83)90226-3)
- [TOCSY paper](https://doi.org/10.1021/ja00295a052)
- [TOCSY paper](https://doi.org/10.1016/0022-2364(85)90018-6)
- [Canonical MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/nmr_liquids/tocsy.m)
- [Spinach Wiki: tocsy.m](https://spindynamics.org/wiki/index.php?title=tocsy.m)
