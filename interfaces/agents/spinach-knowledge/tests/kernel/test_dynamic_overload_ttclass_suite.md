# tests/kernel/test_dynamic_overload_ttclass_suite.m

**Source:** [https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_dynamic_overload_ttclass_suite.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_dynamic_overload_ttclass_suite.m)

## Purpose

Regression test for dynamic dispatch of `ttclass` operator and method overloads. The suite verifies that `ttclass` object operations match exact dense references on one-core and two-core tensor trains, covering constructors, arithmetic, indexing, linear algebra, compression, and solver entry points.

## Behaviour

- Announces the test target with `fprintf('TESTING: Tensor-train overload dispatch\n')` and initialises a result object via `new_test_result('kernel/dynamic_overload_ttclass_suite', ...)`, describing the requirement that `ttclass` overloads match explicit dense references on small deterministic tensor trains.
- Builds one-core tensor trains from `A=[1 2;3 4]` and `B=[0 5;6 1]` as `T=ttclass(2,{A},0)` and `U=ttclass(-3,{B},0)`, with dense references `t_ref=2*A` and `u_ref=-3*B`.
- Checks constructor, `full`, `size`, `sizes`, `ranks`, `numel`, and `subsref` dispatch, including direct function calls `sizes(T)`, `ranks(T)`, and `subsref(T,idx)` with `idx=substruct('()',{2,1})`. Verifies `T.sizes` equals `[2 2]`, `T.ranks` equals `[1;1]`, and structural predicates `T.ncores==1`, `T.ntrains==1`, `ismatrix(T)`, `isnumeric(T)`, `isreal(T)`.
- Checks addition, subtraction, scalar multiplication (left and right `mtimes`), scalar `rdivide`, scalar `mrdivide`, and direct `plus`, `minus`, `rdivide`, `mrdivide` calls against dense references.
- Checks dense-vector and dense-matrix `mtimes` (with `rhs_vec=[1;-2]` and `[1 0;0 -1]`), tensor-train `mtimes`, direct `mtimes`, `dot` (compared to `t_ref'*u_ref`), and `hdot` (compared to the dense Frobenius inner product `sum(conj(t_ref(:)).*u_ref(:))`).
- Builds a complex one-core train `Z=ttclass(1+2i,{[1 2i;-3i 4]},0)` and checks `conj`, transpose (`.'`), `ctranspose` (`'`), direct `ctranspose`, and that `isreal(Z)` is false.
- Checks `trace`, `diag` in both matrix-to-vector and vector-to-matrix directions (using `tt_vec=ttclass(1,{[2;5]},0)`), `sum` and `mean` along dimensions 1 and 2.
- Checks `clearcoeff` (value preserved and coefficient set to 1), `pack` on `T+U`, `ttort` with and without a log-norm output (verifying `exp(lognrm)*full(normalised)` recovers `t_ref` and that `lognrm` is finite), `truncate`, `shrink` on `T+U`, and the Frobenius norm `norm(T,'fro')`.
- Checks direct AMEn summation: `amensum(sum_train,1e-12,sum_opts)` on `ttclass(2,{1;1},0)+ttclass(3,{1;1},0)` with `sum_opts=struct('max_swp',10,'init_guess_rank',1,'enrichment_rank',0,'verb',0)`, expecting the value 5.
- Checks the SPMD-safe save wrapper `save_anyway` by writing `t_ref` to a temporary `.mat` file, loading the `variable` field, deleting the file, and comparing the round-tripped data.
- Checks `unit_like(T)` against `eye(2)`, `rand` on `ttclass(1,{A;B},0)` with rank 2 (seeded via `rng(20240512)`, requiring ranks `[1;2;1]` and finite values), `kron(T,U)` against dense Kronecker multiplication, and `vec(T)` against `t_ref(:)`.
- Builds two-core trains `T2=ttclass(2,{A;B},0)` and `U2=ttclass(-1,{C;D},0)` with `C=[2 -1;0 3]` and `D=[1 0;0 -2]`, checking `full`, two-core `mtimes`, `vec` (against `2*kron(A(:),B(:))`), and `revert` (against `2*kron(B,A)`).
- Checks AMEn solve on a scalar system: `amensolve(lhs,rhs,1e-12,opts,init)` with `lhs=ttclass(2,{1},0)`, `rhs=ttclass(6,{1},0)`, `init=ttclass(1,{1},0)`, and `opts=struct('nswp',2,'verb',0,'enrichment_rank',0,'max_full_size',10)`, expecting the solution 3.
- Verifies that direct `mldivide(T,3)` rejects a non-tensor right-hand side, using `try`/`catch` and `test_true` with a check that the error message contains `'both arguments should be tensor trains'`.
- Numerical comparisons use `test_close` with tolerances of `1e-15` for exact one-core operations, `1e-12` for tensor-train products and compression operations, `1e-10` for two-core multiplication and AMEn solve, and `0` for exact integer/structural comparisons. Structural predicate checks append `PASS` messages directly to `result.messages` and raise errors on failure.

## Inputs and outputs

```matlab
result = test_dynamic_overload_ttclass_suite()
```

- **Outputs:** `result` — regression test result with explanatory messages, as produced by `new_test_result` and accumulated `test_close`/`test_true` calls.
- **Inputs:** None. All test data is deterministic and constructed internally.

## References

- [Spinach MATLAB source (GitHub)](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_dynamic_overload_ttclass_suite.m)
