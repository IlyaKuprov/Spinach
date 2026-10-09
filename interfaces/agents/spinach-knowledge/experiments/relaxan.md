# experiments/relaxan.m

Source: https://github.com/IlyaKuprov/Spinach/blob/main/experiments/relaxan.m
Spinach Wiki: https://spindynamics.org/wiki/index.php?title=relaxan.m

## Purpose

`relaxan` evaluates the configured relaxation model for each spin. It returns longitudinal and transverse rates and times, and the complete relaxation superoperator. This is model analysis, not an experimental measurement.

## Calculation

The routine converts the supplied system to the adjoint representation with `sim2liouv` and calls `relaxation`, with the optional orientation argument when supplied. For each spin it constructs unweighted `Lz` and `L+` vectors with `coil_state`, which remain nonzero even for a zero-population substance, then computes the respective rate as `-real((S'*R*S)/(S'*S))`; each reported time is the reciprocal of its rate. The printed rate columns are labelled Hz and the time columns seconds. Dynamic frequency shifts are dropped.

`euler_angles` is optional and is documented for orientation-dependent relaxation. The source checks for a real three-element value when it is supplied.

## Inputs and outputs

- Call: `[r1,r2,t1,t2,R]=relaxan(spin_system,euler_angles)`; omit the second argument when no orientation is requested.
- `r1` and `r2` are `nspins`-by-1 vectors of longitudinal and transverse rates; `t1` and `t2` are matching time vectors.
- `R` is the complete relaxation superoperator.
- The routine reports each spin number and the system's isotope label alongside the calculated values.

## Source reference

- [Spinach Wiki page](https://spindynamics.org/wiki/index.php?title=relaxan.m)
