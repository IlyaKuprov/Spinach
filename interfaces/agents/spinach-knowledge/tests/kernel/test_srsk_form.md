# tests/kernel/test_srsk_form.m

**Source:** [https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_srsk_form.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_srsk_form.m)

## Purpose

Regression test for the SRSK formalism restriction in Spinach. It verifies that SRSK (source relaxation) requires the spherical-tensor Liouville formalism, and that a Zeeman Liouville request with SRSK is refused explicitly rather than silently selecting a different high-spin relaxation model. The physical scenario is a rapidly relaxing 14N source broadening its scalar-coupled proton.

## Behaviour

The test builds a two-spin system (`1H`, `14N`) at 1e-3 magnet field with a 100 Hz scalar coupling, Lindblad relaxation rates `lind_r1_rates=[1,1e5]` and `lind_r2_rates=[1,1e5]` (rapid nitrogen relaxation), zero equilibrium magnetisation, 298 K temperature, lab-frame relaxation retention, and density-frequency preservation (`rlx_dfs='keep'`).

Under `zeeman-liouv` formalism with `SRSK` in `inter.relaxation` and `inter.srsk_sources=2`, the test loops over all combinations of retention (`labframe`, `secular`, `diagonal`) and equilibrium (`zero`, `IME`, `dibari`) settings, requiring that `relaxation()` throws an error containing `SRSK`, `not implemented`, and `zeeman-liouv`.

It further checks that the refusal remains specific when extended T1/T2 is requested alongside SRSK (`theories={'t1_t2','SRSK'}`), confirming the public SRSK guard precedes the recursive model restriction.

Removing SRSK (`theories={'lindblad'}`) must preserve the specified proton transverse relaxation rate: `R*rho` must equal `-rho` for the `1H` `L+` state within 1e-10 tolerance.

In the supported `sphten-liouv` formalism, the additive SRSK contribution (total relaxation superoperator minus the Lindblad-only superoperator) is checked against Abragam scalar-relaxation rates:

- Longitudinal: `r1_add = (4/3)*coupling^2*(1e-5/(1+freq_diff^2*1e-10))` with `coupling = 2*pi*100`, applied to the `1H` `Lz` state.
- Transverse: `r2_add = (2/3)*coupling^2*(1e-5 + 1e-5/(1+freq_diff^2*1e-10))`, applied to the complex non-Hermitian state `(1+2i)*L+`.

Both checks use 1e-10 absolute and relative tolerances.

## Inputs and outputs

```matlab
result = test_srsk_form()
```

- **Output** `result` — regression test result structure with explanatory messages, created via `new_test_result('kernel/srsk_form', ...)` and accumulated through `test_true` and `test_close` assertions.
- **Input** — none.

## References

- Spinach source: [tests/kernel/test_srsk_form.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_srsk_form.m)
