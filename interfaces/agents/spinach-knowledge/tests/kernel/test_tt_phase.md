# tests/kernel/test_tt_phase.m

## Purpose

Regression test for phase-independent tensor-train compression error budgets. It verifies that TT compression honours absolute error budgets independently of coefficient phase, covering signed and complex coefficients, genuine rank reduction, zero boundaries, and spin products.

Source: [tests/kernel/test_tt_phase.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_tt_phase.m)

## Behaviour

- Initialises a regression record via `new_test_result` with identifier `kernel/tt_phase`, name `Tensor-train coefficient phases`, and the description `Compression must honour absolute error budgets independently of coefficient phase.`
- Saves and restores the MATLAB random number generator state using `onCleanup` for reproducibility; seeds `rng(240924)`.
- Defines the phase set `[1 -1 1i -1i exp(0.37i) 2+3i]` and tracks `max_error` and `reduced_cases`.
- For `ncores` of 2 and 3, builds independent complex cores: each core is `randn(2)+1i*randn(2)` normalised by its Frobenius norm.
- For `nranks` of 1 to 3, uses weights `1`, `[1 1e-5]`, and `[1 0.15 1e-5]` respectively.
- Packs a tolerance-bearing train `P` with per-rank tolerances `1e-3*ones(1,nranks)/nranks` and an exact zero-tolerance train `exact`; computes `reference=full(P)`.
- For scales `0.2` and `7` and each phase, forms `Q=(scale*phase)*P`, applies `shrink`, and compares `full(rounded)` against `expected=(scale*phase)*reference` using `test_close` with tolerance `Q.tolerance` and roundoff `1e-12` under the label `absolute Frobenius truncation budget`.
- For `nranks>1`, asserts the residual exceeds `100*eps*max(1,norm(expected,'fro'))`, asserts all interior ranks of `rounded` are strictly smaller than those of `Q`, and increments `reduced_cases`.
- For zero tolerance, applies `shrink` to `(scale*phase)*exact`, asserts all ranks are retained, and compares with `test_close` using a roundoff of `1e-12*max(1,norm(expected,'fro'))` and zero tolerance under the label `zero tolerance permits only numerical roundoff`.
- Zero-coefficient escape: for tolerances `0` and `1e-3`, shrinks `ttclass(0,cores(:,1),tolerance)` and compares against `zeros(2^ncores)` with zero tolerance and zero roundoff under the label `a zero coefficient must return an exact finite zero`.
- Spin product: defines `spin_x`, `spin_y`, `spin_z` as the standard spin-1/2 operators divided by 2; builds `left=(spin_x+spin_y)/sqrt(2)` and `right=(spin_y+spin_z)/sqrt(2)`; forms `H=ttclass(-2*pi*150,{left;right},1e-8)` and `P=ttclass(1,{eye(2);eye(2)},0)`; expects `-2*pi*150*kron(left,right)`.
- Asserts `norm(expected-expected','fro')==0` (Hermiticity of the expected product), then compares `full(H*P)` with `test_close` at tolerance `1e-8` and roundoff `0` under the label `negative interaction coefficients must survive TT multiplication`, and compares `observed` with `observed'` at tolerance `1e-12` and roundoff `0` under the label `compression must preserve the Hermitian product to roundoff`.
- Prints a summary line: `TT_PHASE max_error=%.15g reduced_cases=%d spin_error=%.15g` with the tracked maximum error, reduced-case count, and `norm(observed-expected,'fro')`.

## Inputs and outputs

- Syntax: `result=test_tt_phase()`
- The function takes no inputs.
- `result` — regression result structure for signed and complex coefficients, genuine rank reduction, zero boundaries, and spin products.

## References

- Source file: [tests/kernel/test_tt_phase.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_tt_phase.m)
