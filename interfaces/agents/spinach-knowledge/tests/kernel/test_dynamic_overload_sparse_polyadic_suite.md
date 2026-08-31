# tests/kernel/test_dynamic_overload_sparse_polyadic_suite.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_dynamic_overload_sparse_polyadic_suite.m`
- Signature: `result=test_dynamic_overload_sparse_polyadic_suite()`
- Total lines: 253

## Purpose

Tests dynamic dispatch of sparse, polyadic, and OPIUM overloads. Syntax: result=test_dynamic_overload_sparse_polyadic_suite()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- The file also defines local helper function(s): `local_have_gpu()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 16-17: Announce the test target; implemented by `fprintf('TESTING: Sparse, polyadic, and OPIUM overload dispatch\n')`.
- Lines 19-22: State the overload target of the test; implemented by `result=new_test_result('kernel/dynamic_overload_sparse_polyadic_suite', 'Dynamic sparse and polyadic overload dispatch', 'RCV, polyadic, and OPIUM overloads must match e…`.
- Lines 24-25: Build deterministic sparse operands; implemented by `A=sparse([1 0 2;0 -3 0])`.
- Lines 30-32: Exercise RCV construction, conversion, and size dispatch; implemented by `result=test_close(result,'rcv sparse conversion',sparse(ra),A,1e-15,1e-15, 'rcv sparse conversion must reproduce the input matrix')`.
- Lines 40-42: Exercise RCV arithmetic dispatch; implemented by `result=test_close(result,'rcv plus rcv',full(ra+rb),full(A+B),1e-15,1e-15, 'rcv plus must add two RCV sparse matrices')`.
- Lines 56-58: Exercise direct RCV method-name dispatch for coverage tracking; implemented by `result=test_close(result,'rcv direct plus',full(plus(ra,rb)),full(A+B),1e-15,1e-15, 'direct plus() must dispatch to the RCV overload')`.
- Lines 68-69: Exercise RCV matrix product and transposition dispatch; implemented by `C=sparse([2 0;-1 3;0 4])`.
- Lines 94-95: Exercise RCV sparsity plotting dispatch offscreen; implemented by `old_visibility=get(groot,'defaultFigureVisible')`.
- Lines 106-107: Exercise guarded RCV GPU dispatch when a GPU is available; implemented by `if local_have_gpu()`.
- Lines 119-121: Build deterministic polyadic operands; implemented by `P=polyadic({{[1 2;0 -1],[0 3;4 5]}, {[2 0;1 1],[-1 2;3 0]}})`.
- Lines 127-128: Exercise polyadic construction, conversion, and structural predicates; implemented by `validate(P)`.
- Lines 139-141: Exercise polyadic arithmetic and Kronecker dispatch; implemented by `result=test_close(result,'polyadic plus',full(P+Q),p_ref+q_ref,1e-15,1e-15, 'polyadic plus must buffer sums without changing values')`.
- Lines 168-170: Exercise polyadic transpose, prefix, suffix, simplify, and inflate dispatch; implemented by `result=test_close(result,'polyadic transpose',full(P.'),p_ref.',1e-15,1e-15, 'polyadic transpose must transpose every represented matrix factor')`.
- Lines 189-190: Exercise guarded polyadic GPU dispatch when a GPU is available; implemented by `if local_have_gpu()`.
- Lines 199-200: Exercise OPIUM construction, conversion, size, multiplication, and kron; implemented by `unit_op=opium(3,1)`.

### Control flow inferred from the code

- Line 107: conditional branch on `local_have_gpu()`.
- Line 109: conditional branch on `~gpu_rcv.isGPU`.
- Line 134: conditional branch on `isempty(P)||~allfinite(P)||nnz(P)~=12`.
- Line 190: conditional branch on `local_have_gpu()`.
- Line 226: conditional branch on `~allfinite(scaled_op)||~isnumeric(scaled_op)||~ismatrix(scaled_op)||nnz(scaled_op)~=1`.

### Key state/data transformations

- Lines 20-22: computes `result` using `result=new_test_result('kernel/dynamic_overload_sparse_polyadic_suite', 'Dynamic sparse and polyadic overload dispatch', 'RCV, polyadic, and OPIUM overloads must match e…`.
- Lines 25: computes `A` using `A=sparse([1 0 2;0 -3 0])`.
- Lines 26: computes `B` using `B=sparse([0 4 -1;5 0 6])`.
- Lines 27: computes `ra` using `ra=rcv(A)`.
- Lines 28: computes `rb` using `rb=rcv(B)`.
- Lines 69: computes `C` using `C=sparse([2 0;-1 3;0 4])`.
- Lines 70: computes `L` using `L=sparse([2 0;-1 3;0 4])`.
- Lines 77: computes `complex_rcv` using `complex_rcv=rcv(A+1i*B)`.
- Lines 95: computes `old_visibility` using `old_visibility=get(groot,'defaultFigureVisible')`.
- Lines 97: computes `plot_cleaner` using `plot_cleaner=onCleanup(@()set(groot,'defaultFigureVisible',old_visibility))`.
- Lines 98: computes `figures_before` using `figures_before=numel(findall(0,'Type','figure'))`.
- Lines 100: computes `figures_after` using `figures_after=numel(findall(0,'Type','figure'))`.
- Lines 108: computes `gpu_rcv` using `gpu_rcv=gpuArray(ra)`.
- Lines 112: computes `cpu_rcv` using `cpu_rcv=gather(gpu_rcv)`.
- Lines 116: computes `result.messages{end+1}` using `result.messages{end+1}='SKIP: rcv gpuArray path skipped because no usable GPU is available.'`.
- Lines 120-121: computes `P` using `P=polyadic({{[1 2;0 -1],[0 3;4 5]}, {[2 0;1 1],[-1 2;3 0]}})`.
- Lines 122: computes `Q` using `Q=polyadic({{[0 1;-2 3],[5 -1;0 2]}})`.
- Lines 123-124: computes `p_ref` using `p_ref=kron([1 2;0 -1],[0 3;4 5])+ kron([2 0;1 1],[-1 2;3 0])`.

### Local helper functions

- Line 234: `local_have_gpu()` — `function answer=local_have_gpu()`. Default to no GPU
  - Representative operation: `answer=false()`.
  - Representative operation: `if exist('gpuDevice','file')~=2`.

## Outputs

- result -regression test result with explanatory messages
- The test exercises object-operation dispatch for rcv, polyadic,
- and opium objects using small deterministic dense references.

## Implementation structure

- Tests dynamic dispatch of sparse, polyadic, and OPIUM overloads. Syntax:
- result=test_dynamic_overload_sparse_polyadic_suite()
- result -regression test result with explanatory messages
- The test exercises object-operation dispatch for rcv, polyadic,
- and opium objects using small deterministic dense references.
- Announce the test target
- State the overload target of the test
- Build deterministic sparse operands
- Exercise RCV construction, conversion, and size dispatch
- Exercise RCV arithmetic dispatch
- Exercise direct RCV method-name dispatch for coverage tracking
- Exercise RCV matrix product and transposition dispatch

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `rcv()`, `test_close()`, `double()`, `plus()`, `minus()`, `times()`, `mtimes()`, `rdivide()`, `ctranspose()`, `horzcat()`, `vertcat()`, `gather()`, `get()`, `set()`, `onCleanup()`.
