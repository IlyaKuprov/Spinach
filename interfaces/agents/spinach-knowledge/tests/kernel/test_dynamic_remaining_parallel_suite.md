# tests/kernel/test_dynamic_remaining_parallel_suite.m

## Purpose

Regression test for the remaining parallel, stochastic, and diagnostic utilities in the kernel layer: distributed-array reconstruction, zero stochastic Redfield integration, and Fokker-Planck overwinding diagnostics.

## Behaviour

- Announces the test target with `TESTING: Remaining parallel and stochastic utilities`.
- Creates a regression test result via `new_test_result` for `kernel/dynamic_remaining_parallel_suite`, described as "Remaining parallel and stochastic utilities", with the requirement that parallel and stochastic helper utilities preserve compact reference cases.
- Ensures a one-worker `Processes` pool exists on the local host via `parpool('Processes',1)` when `gcp('nocreate')` returns empty.
- Distributed-array check: when both `distributed` and `codistributor1d` classes exist (toolbox present), builds `dense_array=reshape(1:12,[3 4])`, calls `distrib_dim(dense_array,2)`, and verifies with `test_close` that `gather(distributed_array)` reproduces the dense array to absolute and relative tolerances of `1e-14`. If the toolbox is unavailable, records a `test_true` pass noting that distributed-array coverage is blocked.
- Redfield integration check: builds a spherical-tensor Liouvillian spin system of dimension 1, sets `H0=sparse(1,1,1e-3,1,1)` and `H1` as 2001 copies of a 1x1 zero sparse matrix, calls `[R,dR]=ngce(spin_system,H0,H1,1,10,0)`, and verifies with `test_close` that both `R` and `dR` equal a 1x1 zero sparse matrix to `1e-14` tolerances (zero relaxation superoperator and zero relaxation uncertainty).
- Overwinding diagnostic check: counts open figures, calls `overwound(rho,[10 1 1],1)` with `rho=ones(10,1)`, closes all figures, and verifies with `test_true` that at least one new figure was drawn.

## Inputs and outputs

```matlab
result=test_dynamic_remaining_parallel_suite()
```

- `result` - regression test result structure with explanatory messages.

No input arguments.

## References

- Source: [tests/kernel/test_dynamic_remaining_parallel_suite.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_dynamic_remaining_parallel_suite.m)
