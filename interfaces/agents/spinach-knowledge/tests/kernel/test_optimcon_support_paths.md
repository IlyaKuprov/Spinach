# tests/kernel/test_optimcon_support_paths.m

**Source:** [https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_optimcon_support_paths.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_optimcon_support_paths.m)

## Purpose

Regression test for small optimal-control support paths in Spinach. It exercises penalty functions, trapezium-product derivatives, objective-function collection, and small line-search helper paths, reporting each check through the standard test-result mechanism.

## Behaviour

- Announces the target with `TESTING: Optimal-control support helper paths` and initialises a result via `new_test_result('optimcon/support_paths', ...)`, describing coverage of `penalty()`, `trapdiff()`, `objeval()`, and line-search helpers.
- Builds a minimal quiet spin system (`sys.output='hush'`, `bas.formalism='zeeman-hilb'`, `tols.small_matrix=64`, `tols.prop_chop=1e-14`) for low-level helper calls.
- **Penalty `none` path:** calls `penalty(waveform,'none',-1,1)` on a 2-by-6 waveform and checks that the value, gradient, and Hessian are all zero.
- **Penalty `NS` path:** checks the norm-square penalty against closed forms — value `sum(waveform(:).^2)/size(waveform,2)`, gradient `2*waveform/size(waveform,2)`, Hessian `2*eye(numel(waveform))/size(waveform,2)` — with tolerances `1e-14`.
- **Penalty `SNS` path:** checks the spillout penalty against explicit clipping residuals `max(waveform-1,0)` and `min(waveform+1,0)`, with value, gradient, and diagonal Hessian references all divided by the number of time points, at tolerance `1e-14`.
- **Penalty `DNS` path:** with bounds `-10,10`, checks that the value is finite, that the gradient matches a centred finite-difference gradient (step `1e-6`, tolerance `1e-6`), and that the Hessian is square over the waveform elements.
- **Penalty `SNSA` path:** on a 4-by-4 Cartesian control waveform with bounds `0,1`, computes channel amplitudes as row-pair Euclidean norms, checks the penalty is positive, matches the value and radial-projection gradient at tolerance `1e-14`, and matches the Hessian against centred finite differences of the gradient at tolerance `1e-6`.
- **`trapdiff` derivatives:** with `dt=1e-3`, `cL=2.0`, `cR=-1.5`, and Pauli-based drift/control Hamiltonians, checks the left- and right-edge trapezium derivative matrices against finite differences of matrix exponentials at tolerance `1e-10`.
- **`objeval` collection:** using a two-channel local objective, checks that value, gradient, and Hessian calls combine channels by subtraction, return consistent objective values, and increment the `fx`, `gfx`, and `hfx` counters appropriately, at tolerance `1e-14`.
- **Line-search conditions:** sets `ls_c1=1e-2`, `ls_c2=0.9`, `ls_tau1=3`, `ls_tau2=0.1`, `ls_tau3=0.5`, then checks `alpha_conds` for monotonic, Armijo, and strong Wolfe curvature acceptance.
- **`cubic_interp`:** checks that the Hermite cubic with opposite endpoint slopes returns maximiser `0.5` and value `0.25` at tolerance `1e-14`.
- **`bracketing`:** on a concave quadratic objective with trial step `0.1`, checks immediate acceptance (`next_act` equal to `'none'`), unchanged step length, objective increase, and gradient shape preservation.
- **`sectioning`:** on the same quadratic, checks that the step length `0.5`, objective `0`, zero gradient, and exit flag `0` are recovered at tolerance `1e-12`.

## Inputs and outputs

**Syntax**

```matlab
result = test_optimcon_support_paths()
```

**Outputs**

- `result` — regression test result with explanatory messages.

## References

- `penalty`, `trapdiff`, `objeval`, `alpha_conds`, `cubic_interp`, `bracketing`, `sectioning` — Spinach optimal-control helper functions exercised by this test.
- `new_test_result`, `test_close`, `test_true` — Spinach kernel test utilities.
- `pauli` — Spinach operator constructor used to build test Hamiltonians.
