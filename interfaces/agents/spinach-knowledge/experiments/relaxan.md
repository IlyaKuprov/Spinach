# experiments/relaxan.m

- Signature: `[r1,r2,t1,t2,R]=relaxan(spin_system,euler_angles)`

## Purpose

Analyses the spin system relaxation model, reports longitudinal and transverse relaxation rates and times for every spin, and returns the relaxation superoperator.

## Numerical / algorithmic content

The routine builds `R` with `relaxation`, optionally at the supplied orientation. For each spin it evaluates the decay rate of its `Lz` state and `L+` state from the corresponding normalized quadratic form with `R`; it returns the rates in Hz and their reciprocals as relaxation times in seconds. Dynamic frequency shifts are dropped.

## Parameters / inputs

- `euler_angles` — optional real three-element Euler-angle vector for orientation-dependent relaxation properties.

## Outputs

- `r1` — vector of longitudinal relaxation rates, in Hz, one per spin.
- `r2` — vector of transverse relaxation rates, in Hz, one per spin.
- `t1` — vector of longitudinal relaxation times, in seconds, one per spin.
- `t2` — vector of transverse relaxation times, in seconds, one per spin.
- `R` — complete relaxation superoperator.

The function requires Liouville-space formalism.
