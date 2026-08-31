# tests/kernel/test_dynamic_rlx_split_suite.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_dynamic_rlx_split_suite.m`
- Signature: `result=test_dynamic_rlx_split_suite()`
- Total lines: 69

## Purpose

Tests relaxation-superoperator component splitting. Syntax: result=test_dynamic_rlx_split_suite()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 16-17: Announce the test target; implemented by `fprintf('TESTING: Relaxation component splitting\n')`.
- Lines 19-22: State the relaxation-splitting target of the test; implemented by `result=new_test_result('kernel/dynamic_rlx_split_suite', 'Relaxation component splitting', 'rlx_split() must partition longitudinal, transverse, and mixed relaxation com…`.
- Lines 24-25: Build a two-spin spherical-tensor Liouville basis; implemented by `sys.magnet=0`.
- Lines 33-34: Interpret the basis using the same documented state categories; implemented by `[L,M]=lin2lm(spin_system.bas.basis)`.
- Lines 40-41: Build a diagonal relaxation matrix with a zero unit-state element; implemented by `matrix_dim=size(spin_system.bas.basis,1)`.
- Lines 46-47: Build the mathematically expected non-overlapping blocks; implemented by `R1_ref=R`.
- Lines 54-55: Split the relaxation superoperator with the production helper; implemented by `[R1,R2,Rm]=rlx_split(spin_system,R)`.
- Lines 57-59: Check the three blocks against their category masks; implemented by `result=test_close(result,'longitudinal block',R1,R1_ref,0,0, 'R1 must contain only single-spin longitudinal states')`.

### Key state/data transformations

- Lines 20-22: computes `result` using `result=new_test_result('kernel/dynamic_rlx_split_suite', 'Relaxation component splitting', 'rlx_split() must partition longitudinal, transverse, and mixed relaxation com…`.
- Lines 25: computes `sys.magnet` using `sys.magnet=0`.
- Lines 26: computes `sys.isotopes` using `sys.isotopes={'1H','1H'}`.
- Lines 27: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={0,0}`.
- Lines 28: computes `inter.temperature` using `inter.temperature=300`.
- Lines 29: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 30: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 31: computes `spin_system` using `spin_system=test_spin_system(sys,inter,bas)`.
- Lines 34: computes `[L,M]` using `[L,M]=lin2lm(spin_system.bas.basis)`.
- Lines 36: computes `mso_mask` using `mso_mask=(sum(logical(spin_system.bas.basis),2)>1)`.
- Lines 41: computes `matrix_dim` using `matrix_dim=size(spin_system.bas.basis,1)`.
- Lines 42: computes `diag_vals` using `diag_vals=(1:matrix_dim)'`.
- Lines 43: computes `diag_vals(~(long_sso_mask|tran_sso_mask|mso_mask))` using `diag_vals(~(long_sso_mask|tran_sso_mask|mso_mask))=0`.
- Lines 44: computes `R` using `R=spdiags(diag_vals,0,matrix_dim,matrix_dim)`.
- Lines 47: computes `R1_ref` using `R1_ref=R`.
- Lines 48: computes `R1_ref(~long_sso_mask,~long_sso_mask)` using `R1_ref(~long_sso_mask,~long_sso_mask)=0`.
- Lines 49: computes `R2_ref` using `R2_ref=R`.
- Lines 50: computes `R2_ref(~tran_sso_mask,~tran_sso_mask)` using `R2_ref(~tran_sso_mask,~tran_sso_mask)=0`.

## Outputs

- result -regression test result with explanatory messages
- The test checks that rlx_split() separates single-spin longitudinal,
- single-spin transverse, and multi-spin relaxation blocks without overlap.

## Implementation structure

- Tests relaxation-superoperator component splitting. Syntax:
- result=test_dynamic_rlx_split_suite()
- result -regression test result with explanatory messages
- The test checks that rlx_split() separates single-spin longitudinal,
- single-spin transverse, and multi-spin relaxation blocks without overlap.
- Announce the test target
- State the relaxation-splitting target of the test
- Build a two-spin spherical-tensor Liouville basis
- Interpret the basis using the same documented state categories
- Build a diagonal relaxation matrix with a zero unit-state element
- Build the mathematically expected non-overlapping blocks
- Split the relaxation superoperator with the production helper

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `rlx_split()`, `test_spin_system()`, `lin2lm()`, `logical()`, `any()`, `diag_vals()`, `spdiags()`, `R1_ref()`, `R2_ref()`, `Rm_ref()`, `test_close()`.
