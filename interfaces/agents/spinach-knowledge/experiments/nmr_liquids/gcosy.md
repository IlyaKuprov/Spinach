# experiments/nmr_liquids/gcosy.m

- Signature: `fid=gcosy(spin_system,parameters,H,R,K)`
- Source: [`experiments/nmr_liquids/gcosy.m`](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/nmr_liquids/gcosy.m)

## Purpose and sequence

This Horne-Morris gradient-selected COSY sequence starts from `Lz` magnetisation on the selected isotope, applies an `Lx` 90-degree pulse, and records the F1 trajectory. Its second pulse has angle `parameters.angle`. A two-gradient sandwich selects the pathway: P uses opposite gradient signs, while N uses equal signs. A positive `parameters.g_stab_del` inserts a stabilisation delay after the first gradient as part of the sandwich propagator and after the second gradient before F2 acquisition. F2 evolution is detected with the selected isotope's `L+` coil state. The source notes P-type selection is less sensitive to mixing-pulse phase errors; P+N returns both pathway components for echo/anti-echo recombination. This is a parameterised simulation sequence, not a measured spectrum or run-verified result.

The Liouvillian is `L=H+1i*R+1i*K`; both dimensions use dwell time `1/parameters.sweep` seconds. Defaults are `g_amp=3` Gauss/cm, `g_dur=2e-3` seconds, `g_stab_del=2e-4` seconds, `s_len=1.5` cm, and `pathway='P'`. They apply when the corresponding field is absent.

## Parameters and inputs

- `parameters.sweep`: positive real scalar sweep width in Hz.
- `parameters.npoints`: two positive integer point counts, ordered F1 then F2.
- `parameters.spins`: one-element cell array naming an isotope present in the system (for example, `{'1H'}` or `{'13C'}`).
- `parameters.angle`: finite real second-pulse angle in radians. The source notes pi/2 is usual and also allows angles such as those for COSY45 and COSY60.
- `parameters.g_amp`: positive real gradient amplitude in Gauss/cm; default 3.
- `parameters.g_dur`: positive real gradient duration in seconds; default 2e-3.
- `parameters.g_stab_del`: non-negative real post-gradient stabilisation delay in seconds; default 2e-4.
- `parameters.s_len`: positive real active sample length in cm; default 1.5.
- `parameters.pathway`: `'P'`, `'N'`, or `'P+N'`; default `'P'`. P uses gradient signs [1,-1], N [1,1].
- `H`, `R`, and `K`: same-sized numeric Hamiltonian, relaxation, and kinetics matrices supplied by the context function. The function requires the `sphten-liouv` formalism.

## Outputs and reference

- `fid`: two-dimensional FID for P or N selection.
- In P+N mode, `fid.pos` is the P-type component and `fid.neg` is the N-type component.
- [Spinach Wiki: `gcosy.m`](https://spindynamics.org/wiki/index.php?title=gcosy.m)
