# tests/kernel/test_dynamic_filter_frontends.m

**Source:** [https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_dynamic_filter_frontends.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_dynamic_filter_frontends.m)

## Purpose

Regression test for the dynamic state-filter front-end kernels of Spinach: `coherence()`, `correlation()`, `decouple()`, `homospoil()`, and `spinlock()`. The test runs these functions on small spherical-tensor Liouville-space spin systems whose surviving state components are analytically known, and checks that each filter keeps exactly the intended basis components.

## Behaviour

- Registers a test target named `kernel/dynamic_filter_frontends` ("Dynamic state-filter front ends") via `new_test_result()`, with the requirement that state-selection front ends keep exactly the intended basis components.
- **Coherence and correlation:** builds a heteronuclear `1H`/`13C` spherical-tensor Liouville system (`sys.magnet=14.1`, Zeeman scalars `{1.0,2.0}`, scalar coupling `10.0` between spins 1 and 2, formalism `'sphten-liouv'`, approximation `'none'`) through `test_spin_system()`. It constructs single- and two-spin state components (`Lz` on each spin, `L+`/`L-` on the proton combined with `Lz` on carbon, and `Lz`-`Lz`), then:
  - Applies `coherence()` with selection `{{'1H',1}}` to a two-column stack of mixed plus/minus proton coherences and checks that only the proton plus-one coherence order is retained.
  - Applies `correlation()` with order `1` and `'all'` to a mix of one-spin and two-spin terms and checks that only the one-spin terms survive.
  - Applies `correlation()` with order `1` and numeric spin index `2` and checks that only carbon-local terms survive (the numeric spin-list path).
- **Decoupling:** builds the same heteronuclear system, applies `assume(spin_system,'nmr')`, and computes the Hamiltonian with `hamiltonian()`. From a state containing proton, carbon, and proton-carbon terms, it calls `decouple()` with `{'13C'}` and with numeric index `2`, checking in both cases that only the proton `Lz` component remains. It also verifies that the rows and columns of the returned Hamiltonian touching carbon (basis column 2 nonzero) are zeroed, using 1-norms of the masked rows and columns.
- **Homospoil:** builds a homonuclear two-`1H` system with zero Zeeman scalars and zero couplings, mixes longitudinal (`Lz`-`Lz`), zero-quantum (`L+`-`L-`), and double-quantum (`L+`-`L+`) components, then checks that `homospoil(...,'keep')` retains the longitudinal and zero-quantum terms while `homospoil(...,'destroy')` retains only the longitudinal spherical tensors.
- **Spinlock:** builds a one-spin `1H` system with zero Zeeman scalar, applies `assume(spin_system,'nmr')`, prepares `Lx`, `Ly` operators and a mixed state `rho_x + 2*rho_y + 3*rho_z`, then checks that `spinlock(...,'X')` returns only the X component and `spinlock(...,'Y')` returns only `2*rho_y`.
- All comparisons use `test_close()`; coherence, correlation, decoupling, and homospoil checks use tolerances `1e-14` (absolute and relative), while spinlock checks use `1e-13`.

## Inputs and outputs

```matlab
result = test_dynamic_filter_frontends()
```

- **Output:** `result` — regression test result structure with explanatory messages, accumulated by the `test_close()` checks.
- **Input:** none.

## References

- Tested functions: `coherence()`, `correlation()`, `decouple()`, `homospoil()`, `spinlock()`.
- Supporting test infrastructure: `new_test_result()`, `test_close()`, `test_spin_system()`, `state()`, `operator()`, `hamiltonian()`, `assume()`.
