# tests/kernel/test_matrix_utility_suite.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_matrix_utility_suite.m`
- Signature: `result=test_matrix_utility_suite()`
- Total lines: 118

## Purpose

Tests small matrix utility functions. Syntax: result=test_matrix_utility_suite()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 16-17: Announce the test target; implemented by `fprintf('TESTING: Matrix utility functions\n')`.
- Lines 19-22: State the utility target of the test; implemented by `result=new_test_result('kernel/matrix_utility_suite', 'Matrix utility functions', 'matrix helpers must preserve their exact algebraic definitions.')`.
- Lines 24-25: Define small test matrices; implemented by `A=[1 2;3 4]`.
- Lines 28-30: Check anticommutator and cheap norm; implemented by `result=test_close(result,'acomm',acomm(A,B),A*B+B*A,1e-15,1e-15, 'the anticommutator is AB+BA')`.
- Lines 34-35: Check polyadic norm estimation paths; implemented by `P=polyadic({{A}})`.
- Lines 41-42: Check rectangular polyadic sign-history dimensions; implemented by `R=[1 0;2 3;4 5]`.
- Lines 51-52: Check the Higham-Tisseur zero-sign convention; implemented by `Z=[-2 0;-1 2]`.
- Lines 57-58: Check complex adjoint and phase handling; implemented by `Z=[-2 0;-1+1i 2]`.
- Lines 63-64: Check the block estimator lower-bound contract; implemented by `rng_state=rng`.
- Lines 74-76: Check polyadic realness dispatch; implemented by `result=test_true(result,'polyadic isreal true',isreal(polyadic({{Z}})), 'real polyadic cores are recognised without opening the Kronecker products')`.
- Lines 80-82: Check identity and trace predicates; implemented by `result=test_true(result,'iseye true',iseye(speye(3)), 'the sparse unit matrix is recognised as identity')`.
- Lines 94-95: Check sparse block-diagonal assembly; implemented by `S=sp_block_diag([1 2;3 4],[5;6])`.
- Lines 100-101: Check row and column replication; implemented by `M=[1 2 3;4 5 6]`.
- Lines 107-108: Check sparse cleanup drops sub-tolerance elements; implemented by `spin_system.sys.output='hush'`.

### Key state/data transformations

- Lines 20-22: computes `result` using `result=new_test_result('kernel/matrix_utility_suite', 'Matrix utility functions', 'matrix helpers must preserve their exact algebraic definitions.')`.
- Lines 25: computes `A` using `A=[1 2;3 4]`.
- Lines 26: computes `B` using `B=[0 1;-1 2]`.
- Lines 35: computes `P` using `P=polyadic({{A}})`.
- Lines 42: computes `R` using `R=[1 0;2 3;4 5]`.
- Lines 52: computes `Z` using `Z=[-2 0;-1 2]`.
- Lines 64: computes `rng_state` using `rng_state=rng`.
- Lines 65: computes `rng_cleanup` using `rng_cleanup=onCleanup(@()rng(rng_state))`.
- Lines 69: computes `est` using `est=cheap_norm(P,3,5)`.
- Lines 95: computes `S` using `S=sp_block_diag([1 2;3 4],[5;6])`.
- Lines 96: computes `S_ref` using `S_ref=sparse([1 2 0;3 4 0;0 0 5;0 0 6])`.
- Lines 101: computes `M` using `M=[1 2 3;4 5 6]`.
- Lines 108: computes `spin_system.sys.output` using `spin_system.sys.output='hush'`.
- Lines 109: computes `spin_system.sys.disable` using `spin_system.sys.disable={}`.
- Lines 110: computes `spin_system.sys.enable` using `spin_system.sys.enable={}`.
- Lines 111: computes `spin_system.tols.dense_matrix` using `spin_system.tols.dense_matrix=0.9`.
- Lines 112: computes `spin_system.tols.small_matrix` using `spin_system.tols.small_matrix=10`.
- Lines 113: computes `C` using `C=sparse([1e-12 1;0 2])`.

## Outputs

- result -regression test result with explanatory messages
- The test checks low-level matrix helpers that Spinach uses throughout the
- kernel for algebra, sparsity, block assembly, and indexing operations.

## Implementation structure

- Tests small matrix utility functions. Syntax:
- result=test_matrix_utility_suite()
- result -regression test result with explanatory messages
- The test checks low-level matrix helpers that Spinach uses throughout the
- kernel for algebra, sparsity, block assembly, and indexing operations.
- Announce the test target
- State the utility target of the test
- Define small test matrices
- Check anticommutator and cheap norm
- Check polyadic norm estimation paths
- Check rectangular polyadic sign-history dimensions
- Check the Higham-Tisseur zero-sign convention

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `test_close()`, `acomm()`, `cheap_norm()`, `polyadic()`, `onCleanup()`, `rng()`, `test_true()`, `clear()`, `iseye()`, `speye()`, `istraceless()`, `krondelta()`, `sp_block_diag()`, `repcols()`, `reprows()`.
