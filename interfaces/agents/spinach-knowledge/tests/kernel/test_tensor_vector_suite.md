# tests/kernel/test_tensor_vector_suite.m

## Purpose

Regression test suite for Spinach kernel tensor, vector, distribution, and relaxation utilities. The suite verifies Hermite spline interpolation, the skew-normal density in its normal limit, Fokker-Planck vector reshaping helpers, rotational correlation-function coefficients, tensor isotope-shift helpers, and small spin-system tensor extractors against exact analytical cases.

## Numerical invariants

- **Interpolation and probability:** `herm_spline` exactly reproduces a quadratic from its endpoint values and slopes (`1e-15` tolerance). At zero skew, `snormpdf` reduces to the corresponding normal density (`1e-14`).
- **Spatial Fokker–Planck embedding:** a small phantom and spin state combine as `kron(phantom(:),spin_state)`. Projection onto the coil recovers the phantom, while spatial averaging returns `mean(phantom(:))*spin_state`; all three identities are checked at `1e-15`.
- **Rotational correlation:** for isotropic correlation time `tau_c=2e-9` s, the first `corrfun` weight is `1/5`, its rate is `-1/tau_c`, and its species state mask is `[1;1;1;0]` (rate tolerance `1e-12`).
- **Spin tensors and offsets:** replacing the isotropic shielding changes `trace(tensor)/3` to `5` without changing its anisotropy (`1e-13`). `get_coupling` adds both directed matrices into `diag([5 7 9])`; `gtensorof` follows `-ddscal*gamma*hbar/muB` and gives `-4*eye(3)` in the reference case. `offsetof` converts the residual angular frequency by `-1/(2*pi)` to `-5`. The remaining matrix comparisons use `1e-14`–`1e-15`.

## Inputs and outputs

**Syntax**

```matlab
result=test_tensor_vector_suite()
```

**Outputs**

- `result` — regression test result object with explanatory messages accumulated across all checks.

The function takes no inputs.

## References

- Source: [tests/kernel/test_tensor_vector_suite.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_tensor_vector_suite.m)
