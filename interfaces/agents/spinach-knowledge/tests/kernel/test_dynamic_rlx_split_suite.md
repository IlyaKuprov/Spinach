# tests/kernel/test_dynamic_rlx_split_suite.m

Source: [https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_dynamic_rlx_split_suite.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_dynamic_rlx_split_suite.m)

## Purpose

Regression test for `rlx_split()`, the relaxation-superoperator component splitter. It verifies that the function partitions a relaxation superoperator into single-spin longitudinal, single-spin transverse, and multi-spin (mixed) relaxation blocks without overlap.

## Behaviour

- Announces the target with `fprintf('TESTING: Relaxation component splitting\n')` and initialises a test result via `new_test_result()` under the identifier `kernel/dynamic_rlx_split_suite`.
- Builds a two-spin `1H` spherical-tensor Liouville system (`bas.formalism='sphten-liouv'`, `bas.approximation='none'`, `inter.temperature=300`, zero scalar Zeeman couplings) using `test_spin_system()`.
- Derives state-category masks from the basis with `lin2lm()`:
  - `sso_mask`: single-spin states (`sum(logical(basis),2)==1`).
  - `mso_mask`: multi-spin states (`sum(logical(basis),2)>1`).
  - `long_sso_mask`: single-spin longitudinal states (`L>0`, `M==0`).
  - `tran_sso_mask`: single-spin transverse states (`L>0`, `M~=0`).
- Constructs a diagonal relaxation matrix `R` of dimension `size(spin_system.bas.basis,1)` with diagonal values `1:matrix_dim`, zeroed for elements outside all three category masks.
- Builds reference blocks `R1_ref`, `R2_ref`, `Rm_ref` by masking `R` to the longitudinal, transverse, and multi-spin categories respectively.
- Calls `[R1,R2,Rm]=rlx_split(spin_system,R)` and compares each block to its reference with `test_close()` using tolerances `0,0`:
  - `R1` must contain only single-spin longitudinal states.
  - `R2` must contain only single-spin transverse states.
  - `Rm` must contain only multi-spin, or otherwise mixed, states.
  - `R1+R2+Rm` must exactly reconstruct the diagonal category-partitioned matrix `R`.

## Inputs and outputs

```matlab
result=test_dynamic_rlx_split_suite()
```

- **Outputs**: `result` — regression test result structure with explanatory messages.
- **Inputs**: none.

## References

- `rlx_split` — relaxation-superoperator component splitter under test.
- `new_test_result`, `test_close` — test harness utilities.
- `test_spin_system` — test spin-system constructor.
- `lin2lm` — converts basis elements to spherical-tensor indices `L` and `M`.
