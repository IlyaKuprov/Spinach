# tests/kernel/test_remtrace_tensor.m

**Source:** https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_remtrace_tensor.m

## Purpose

Regression test for the `remtrace` kernel function. It verifies that `remtrace` subtracts the isotropic component of a second-rank interaction tensor, leaving the anisotropic traceless part.

## Behaviour

- Announces the test target with `fprintf('TESTING: Traceless rank-two tensor construction\n')`.
- Registers a new test result via `new_test_result` with suite name `kernel/remtrace_tensor`, description `'Traceless rank-two tensor construction'`, and the physical statement `'anisotropic interaction tensors are obtained by subtracting the isotropic trace.'`.
- Defines a symmetric interaction tensor with non-zero isotropic part:

  ```matlab
  A = [1 2 0; 2 3 0; 0 0 5];
  ```

- Computes the reference traceless tensor as `A_ref = A - eye(3)*trace(A)/3` and the observed result as `A_obs = remtrace(A)`.
- Runs three closeness checks with tolerances `1e-15` (absolute and relative):
  1. `'explicit isotropic subtraction'` — compares `A_obs` against `A_ref`, with the message `'the isotropic part is trace(A)/3 times the unit matrix'`.
  2. `'zero trace'` — compares `trace(A_obs)` against `0`, with the message `'anisotropic second-rank tensors have zero trace'`.
  3. `'anisotropy preserved'` — compares `A_obs(1,2)` against `A(1,2)`, with the message `'subtracting isotropic trace leaves off-diagonal anisotropy unchanged'`.

## Inputs and outputs

```matlab
result = test_remtrace_tensor()
```

- **Outputs:**
  - `result` — regression test result structure with explanatory messages, accumulated by `new_test_result` and `test_close`.
- **Inputs:** none.

## References

- `remtrace` — kernel function under test (trace removal from a second-rank tensor).
- `new_test_result`, `test_close` — test harness helpers used to register and evaluate checks.
