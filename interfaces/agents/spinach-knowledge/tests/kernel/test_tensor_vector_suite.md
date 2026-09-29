# tests/kernel/test_tensor_vector_suite.m

## Purpose

Regression test suite for Spinach kernel tensor, vector, distribution, and relaxation utilities. The suite verifies Hermite spline interpolation, the skew-normal density in its normal limit, Fokker-Planck vector reshaping helpers, rotational correlation-function coefficients, tensor isotope-shift helpers, and small spin-system tensor extractors against exact analytical cases.

## Behaviour

The function announces the test target with `fprintf`, initialises a regression test result via `new_test_result` for `kernel/tensor_vector_suite`, and then runs a sequence of checks, each appending pass/fail information with explanatory messages:

- **`herm_spline` quadratic reproduction**: interpolates on `linspace(0,1,6)` with endpoint values `0` and `1` and derivatives `0` and `2`, comparing against `spline_grid.^2` with absolute and relative tolerances `1e-15`.
- **`snormpdf` zero-skew normal limit**: evaluates on the grid `-2:2` with location `1`, scale `2`, and skewness `0`, comparing against the ordinary normal density `exp(-0.5*((pdf_grid-1)/2).^2)/(2*sqrt(2*pi))` with tolerances `1e-14`.
- **`phan2fpl` Kronecker embedding**: embeds the phantom `[1 3;2 4]` with spin state `[1;2]`, comparing against `kron(phantom(:),spin_state)` with tolerances `1e-15`.
- **`fpl2phan` observable projection**: extracts the painted image from the Fokker-Planck vector using coil `[1;0]` and grid dimensions `[2 2]`, recovering the original phantom with tolerances `1e-15`.
- **`fpl2rho` spatial average**: averages the Fokker-Planck vector over spatial cells `[2 2]`, comparing against `mean(phantom(:))*spin_state` with tolerances `1e-15`.
- **`corrfun` isotropic coefficients**: uses a minimal `sphten-liouv` spin system with basis `[1 0;0 1;1 1;0 0]`, chemical species partition `{1:2}`, and `rlx.tau_c={2e-9}`; calls `corrfun(corr_system,2,3,3,3,3)` and checks that the first weight equals `1/5` (tolerance `1e-15`), the first rate equals `-1/tau_c` (tolerance `1e-12`), and the first state mask equals `logical([1;1;1;0])`.
- **`shift_iso` isotropic replacement**: applies `shift_iso({diag([1 2 6])},1,5)` and checks that the trace divided by 3 equals `5` (tolerance `1e-13`) and that the anisotropic part `shifted - eye(3)*trace(shifted)/3` equals the original anisotropy `orig_tensor - eye(3)*trace(orig_tensor)/3` (tolerance `1e-13`).
- **`get_coupling` bidirectional sum**: populates `inter.coupling.matrix{1,2}=diag([1 2 3])` and `inter.coupling.matrix{2,1}=diag([4 5 6])` on a two-spin system (`1H`, `13C`), checking that `get_coupling(spin_system,1,2)` equals `diag([5 7 9])` with tolerances `1e-15`.
- **`gtensorof` scaling conversion**: sets `inter.zeeman.ddscal={4*eye(3)}`, `inter.gammas=2`, `tols.hbar=3`, `tols.muB=6`, checking that `gtensorof(spin_system,1)` equals `-4*eye(3)` with tolerances `1e-15`, consistent with the relation `-ddscal*gamma*hbar/muB`.
- **`offsetof` isotropic shift**: sets `inter.zeeman.matrix={2*pi*20*eye(3)}` and `inter.basefrqs=2*pi*15`, checking that `offsetof(spin_system,1)` equals `-5` with tolerances `1e-14`, consistent with the negative residual angular frequency divided by `2*pi`.

Two local helper functions construct minimal spin systems: `local_corr_system` builds the `sphten-liouv` system used by `corrfun`, and `local_tensor_system` builds a two-spin `1H`/`13C` system with empty `inter` and `tols` structs for the tensor extractor tests.

## Inputs and outputs

**Syntax**

```matlab
result = test_tensor_vector_suite()
```

**Outputs**

- `result` — regression test result object with explanatory messages accumulated across all checks.

The function takes no inputs.

## References

- Source: [tests/kernel/test_tensor_vector_suite.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_tensor_vector_suite.m)
