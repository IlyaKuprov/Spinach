# tests/kernel/test_dynamic_overload_ttclass_suite.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_dynamic_overload_ttclass_suite.m`
- Signature: `result=test_dynamic_overload_ttclass_suite()`
- Total lines: 216

## Purpose

Tests dynamic dispatch of tensor-train class overloads. Syntax: result=test_dynamic_overload_ttclass_suite()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 16-17: Announce the test target; implemented by `fprintf('TESTING: Tensor-train overload dispatch\n')`.
- Lines 19-22: State the overload target of the test; implemented by `result=new_test_result('kernel/dynamic_overload_ttclass_suite', 'Dynamic ttclass overload dispatch', 'ttclass overloads must match explicit dense references on small det…`.
- Lines 24-25: Build one-core tensor trains and dense references; implemented by `A=[1 2;3 4]`.
- Lines 32-34: Exercise constructor, full, size, sizes, ranks, numel, and subsref dispatch; implemented by `result=test_close(result,'ttclass full one-core',full(T),t_ref,1e-15,1e-15, 'ttclass full must open a one-core tensor train')`.
- Lines 59-61: Exercise addition, subtraction, scalar multiplication, and scalar division; implemented by `result=test_close(result,'ttclass plus',full(T+U),t_ref+u_ref,1e-15,1e-15, 'ttclass plus must concatenate buffered trains without changing values')`.
- Lines 81-82: Exercise matrix, vector, and tensor-train multiplication dispatch; implemented by `rhs_vec=[1;-2]`.
- Lines 96-97: Exercise conjugation, transposition, trace, diagonal, sum, and mean dispatch; implemented by `Z=ttclass(1+2i,{[1 2i;-3i 4]},0)`.
- Lines 127-128: Exercise coefficient clearing, packing, orthogonalisation, truncation, and shrinkage; implemented by `cleared=clearcoeff(T)`.
- Lines 156-157: Exercise direct AMEn summation and isolated SPMD-save wrapper dispatch; implemented by `sum_opts=struct('max_swp',10,'init_guess_rank',1,'enrichment_rank',0,'verb',0)`.
- Lines 169-171: Exercise unit_like, rand, kron, and vec on object instances; implemented by `result=test_close(result,'ttclass unit_like',full(unit_like(T)),eye(2),1e-15,1e-15, 'unit_like must build an identity tensor train with matching topology')`.
- Lines 183-184: Exercise two-core full, multiplication, vectorisation, and bit-reversal paths; implemented by `C=[2 -1;0 3]`.
- Lines 199-200: Exercise AMEn solve on the smallest deterministic tensor train system; implemented by `lhs=ttclass(2,{1},0)`.

### Control flow inferred from the code

- Line 54: conditional branch on `(T.ncores~=1)||(T.ntrains~=1)||~ismatrix(T)||~isnumeric(T)||~isreal(T)`.
- Line 107: conditional branch on `isreal(Z)`.
- Line 143: conditional branch on `~isfinite(lognrm)`.
- Line 174: conditional branch on `~isa(rand_train,'ttclass')||any(rand_train.ranks~=[1;2;1])||~all(isfinite(full(rand_train)),'all')`.

### Key state/data transformations

- Lines 20-22: computes `result` using `result=new_test_result('kernel/dynamic_overload_ttclass_suite', 'Dynamic ttclass overload dispatch', 'ttclass overloads must match explicit dense references on small det…`.
- Lines 25: computes `A` using `A=[1 2;3 4]`.
- Lines 26: computes `B` using `B=[0 5;6 1]`.
- Lines 27: computes `T` using `T=ttclass(2,{A},0)`.
- Lines 28: computes `U` using `U=ttclass(-3,{B},0)`.
- Lines 29: computes `t_ref` using `t_ref=2*A`.
- Lines 30: computes `u_ref` using `u_ref=-3*B`.
- Lines 51: computes `idx` using `idx=substruct('()',{2,1})`.
- Lines 57: computes `result.messages{end+1}` using `result.messages{end+1}='PASS: ttclass ncores, ntrains, ismatrix, isnumeric, and isreal predicates returned expected values.'`.
- Lines 82: computes `rhs_vec` using `rhs_vec=[1;-2]`.
- Lines 97: computes `Z` using `Z=ttclass(1+2i,{[1 2i;-3i 4]},0)`.
- Lines 98: computes `z_ref` using `z_ref=full(Z)`.
- Lines 115: computes `tt_vec` using `tt_vec=ttclass(1,{[2;5]},0)`.
- Lines 128: computes `cleared` using `cleared=clearcoeff(T)`.
- Lines 133: computes `packed` using `packed=pack(T+U)`.
- Lines 136: computes `ortho` using `ortho=ttort(T,+1)`.
- Lines 139: computes `[normalised,lognrm]` using `[normalised,lognrm]=ttort(T,+1)`.
- Lines 147: computes `truncated` using `truncated=truncate(ortho)`.

## Outputs

- result -regression test result with explanatory messages
- The test exercises ttclass object operations against exact dense
- references on one-core and two-core tensor trains.

## Implementation structure

- Tests dynamic dispatch of tensor-train class overloads. Syntax:
- result=test_dynamic_overload_ttclass_suite()
- result -regression test result with explanatory messages
- The test exercises ttclass object operations against exact dense
- references on one-core and two-core tensor trains.
- Announce the test target
- State the overload target of the test
- Build one-core tensor trains and dense references
- Exercise constructor, full, size, sizes, ranks, numel, and subsref dispatch
- Exercise addition, subtraction, scalar multiplication, and scalar division
- Exercise matrix, vector, and tensor-train multiplication dispatch
- Exercise conjugation, transposition, trace, diagonal, sum, and mean dispatch

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `ttclass()`, `test_close()`, `double()`, `sizes()`, `ranks()`, `t_ref()`, `substruct()`, `subsref()`, `ismatrix()`, `plus()`, `minus()`, `rdivide()`, `mrdivide()`, `mtimes()`, `dot()`.
