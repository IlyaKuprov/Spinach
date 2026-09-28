# tests/kernel/test_dynamic_filter_frontends.m

- Signature: `result=test_dynamic_filter_frontends()`

## Purpose

Tests dynamic state-filter front-end kernels on compact spin systems with analytically known surviving state components.

## Physical / mathematical content

The checks use spherical-tensor Liouville-space systems: a coupled `1H`–`13C` pair for coherence, correlation, and decoupling filters; a two-`1H` system for longitudinal, zero-quantum, and double-quantum components; and a one-`1H` system for transverse spin-lock projections.

## Numerical / algorithmic content

- `coherence()` retains the requested proton `+1` coherence from mixed `L+`/`L-` state columns.
- `correlation(...,1,'all')` retains one-spin terms and removes the two-spin term. With numeric spin list `2`, it retains the carbon-local and proton–carbon terms and removes the proton-local term.
- `decouple()` removes state components involving carbon when carbon is specified by either isotope name or numeric index. The isotope-name check also verifies that Liouvillian rows and columns whose basis components involve carbon are zero.
- `homospoil(...,'keep')` retains longitudinal and zero-quantum terms but removes the double-quantum term; `'destroy'` retains only the longitudinal term.
- `spinlock(...,'X')` retains the X component and removes Y and Z; `'Y'` retains the Y component and removes X and Z.

## Outputs

- `result` — regression test result with explanatory messages. The comparisons use absolute and relative tolerances of `1e-14` for the coherence, correlation, decoupling, and homospoil checks, and `1e-13` for the spin-lock checks.

## Attribution

- ilya.kuprov@weizmann.ac.il