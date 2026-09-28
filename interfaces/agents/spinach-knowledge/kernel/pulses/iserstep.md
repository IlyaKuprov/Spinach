# kernel/pulses/iserstep.m

- Signature: `rho_b=iserstep(spin_system,LTM,rho_a,dt)`

## Purpose

Lie-group and Runge–Kutta–Munthe-Kaas solvers for the Lie equation. The LG methods implement Equation A.1, with minor typos corrected, from [the cited paper](http://dx.doi.org/10.1088/0305-4470/39/19/S07). Unlike `step()`, the Liouvillian may depend on the density matrix.



## Parameters / inputs

- `spin_system` — Spinach data structure from the `create.m` and `basis.m` constructors.
- `L` — function handle `L(t,rho)` that takes time and state vector and returns the evolution generator (rad/s) for `d_rho/d_t = -i*L(t,rho)*rho`.
- `rho_a` — state vector at the start of the evolution period.
- `t` — time at the start of the evolution, in seconds.
- `dt` — evolution time step, in seconds.
- `method` — one of `'PWCL'`, `'LG2'`, `'LG4'`, `'LG4A'`, `'RKMK4'`, `'RKMK-DP5'`, `'RKMK-DP8'`, or `'RKMK-RKF45'`; `'LG4'` has a good balance of efficiency and numerical accuracy.

## Outputs

- `rho_b` — state vector at the end of the evolution time step.
