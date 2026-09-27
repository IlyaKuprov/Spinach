# experiments/nmr_liquids/gcosy.m

- Signature: `fid=gcosy(spin_system,parameters,H,R,K)`

## Purpose

Horne-Morris gradient-selected COSY pulse sequence.

## Implementation

The sequence applies an initial 90-degree `Lx` pulse and evolves the F1 trajectory. It forms the second-pulse propagator using `parameters.angle`, then uses a gradient sandwich for pathway selection: `P` uses opposite gradient signs, while `N` uses equal signs. Optional `parameters.g_stab_del` delays are applied after the first and second gradients. In `P+N` mode both pathways are returned in `fid.pos` and `fid.neg` for echo/anti-echo recombination; otherwise the selected pathway is returned as a two-dimensional `fid`. The documented default `P` pathway is less sensitive to mixing-pulse phase errors.

## Parameters / inputs

- `parameters.sweep` — sweep width in Hz.
- `parameters.npoints` — number of points for both dimensions.
- `parameters.spins` — nuclei on which the sequence runs, specified as `{'1H'}`, `{'13C'}`, etc.
- `parameters.angle` — second pulse angle in radians, usually `pi/2`, but also allows COSY45, COSY60, etc.
- `parameters.g_amp` — gradient amplitude in Gauss/cm; defaults to 3.
- `parameters.g_dur` — gradient duration in seconds; defaults to `2e-3`.
- `parameters.g_stab_del` — post-gradient stabilisation delay in seconds; defaults to `2e-4`.
- `parameters.s_len` — active sample length in cm; defaults to 1.5.
- `parameters.pathway` — optional coherence pathway selection, either `'P'`, `'N'`, or `'P+N'`; defaults to `'P'`.
- `H` — Hamiltonian matrix, received from context function.
- `R` — relaxation superoperator, received from context function.
- `K` — kinetics superoperator, received from context function.

The implementation requires `sphten-liouv` formalism, same-sized matrix inputs `H`, `R`, and `K`, and two positive integer point counts.

## Outputs

- `fid` — two-dimensional free induction decay, or a structure with P-type `fid.pos` and N-type `fid.neg` fields in `P+N` mode.

## References

- [Spin Dynamics Wiki: `gcosy.m`](https://spindynamics.org/wiki/index.php?title=gcosy.m)
