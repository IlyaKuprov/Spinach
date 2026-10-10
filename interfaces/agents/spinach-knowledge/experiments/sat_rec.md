# experiments/sat_rec.m

Source: https://github.com/IlyaKuprov/Spinach/blob/main/experiments/sat_rec.m
Spinach Wiki: https://spindynamics.org/wiki/index.php?title=sat_rec.m

## Purpose and initial condition

`sat_rec` computes a saturation-recovery sequence with analytical saturation: it sets the initial state to `unit_state(spin_system)` rather than applying an explicit saturation pulse. The source notes that the relaxation superoperator must be thermalised. It does not accept a user-supplied starting state.

## Propagation and acquisition

The routine forms `L=H+1i*R+1i*K`. It generates a relaxation trajectory from the unit state using a step of `max_delay/n_delays` for `n_delays` trajectory steps. At every trajectory state it applies a `pi/2` pulse about the Y component of `L+` for the selected isotope, `parameters.spins{1}`. It then acquires using the corresponding `L+` detection state and the same propagation generator, with dwell `1/sweep` and `npoints-1` evolution steps.

The `sweep` parameter is documented in Hz, so the FID dwell is `1/sweep`. `max_delay` is documented as the longest relaxation delay; the source does not state a unit in its parameter description. The relaxation trajectory includes its initial zero-delay state, so observable evolution returns `fids` with `npoints` time-sample rows and `n_delays+1` columns, one FID per delay including zero.

## Required parameters and outputs

- `sweep`: positive real scalar, Hz; `npoints`: positive integer.
- `spins`: one-element cell array containing an isotope present in the system, e.g. `{'1H'}` (or `{'13C'}`).
- `max_delay`: positive real scalar, the longest recovery delay; `n_delays`: positive integer.
- `H`, `R` and `K`: numeric matrices with matching dimensions. The source requires the relaxation superoperator to be thermalised and rejects a system whose `spin_system.rlx.equilibrium` is `'zero'`.
- `fids`: `npoints × (n_delays+1)` matrix, each column a FID for one relaxation delay including zero.

This describes the implemented sequence, not simulated or measured recovery curves.

## Source reference

- [Spinach Wiki page](https://spindynamics.org/wiki/index.php?title=sat_rec.m)
