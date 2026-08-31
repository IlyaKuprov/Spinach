# examples/fundamentals/tensor_structures/polyadic_test_2.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/fundamentals/tensor_structures/polyadic_test_2.m`
- Signature: `polyadic_test_2()`
- Total lines: 110

## Purpose

Unit tests for advanced polyadic functionality.

## Physical / mathematical content

- Fundamentals examples. These are unit tests, convention checks, and pedagogical demonstrations of operator algebra, perturbation theory, tensor conventions, symmetry, quadrature, and numerical differentiation.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- The file also defines local helper function(s): `assert_small()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 7-8: Get random test matrices; implemented by `a=randn(2,2)+1i*randn(2,2)`.
- Lines 17-18: Build a reference polyadic and its matrix form; implemented by `p=polyadic({{a,b},{c,d}})`.
- Lines 21-22: Check constructor, full, inflate, and validate; implemented by `assert_small(norm(full(p)-full(ref),1),1e-12)`.
- Lines 27-28: Check prefixes, suffixes, size, and emptiness; implemented by `p_pref=prefix(left,p)`.
- Lines 36-37: Check addition and subtraction paths; implemented by `q=polyadic({{d.',a.'}})`.
- Lines 44-45: Check multiplication paths; implemented by `assert_small(norm(full(3*p)-full(3*ref),1),1e-12)`.
- Lines 56-57: Check Kronecker products; implemented by `k_num=randn(2,2)+1i*randn(2,2)`.
- Lines 63-64: Check transpose operations; implemented by `assert_small(norm(full(transpose(p))-full(transpose(ref)),1),1e-12)`.
- Lines 68-69: Check finiteness and internal non-zero counts; implemented by `assert(isequal(allfinite(p),all(isfinite(ref(:)))))`.
- Lines 75-76: Check zero-dimension behaviour; implemented by `p_empty=polyadic({{zeros(0,2)}})`.
- Lines 81-82: Check nested simplification paths; implemented by `p_nested=polyadic({{polyadic({{eye(2)}}),b}})`.
- Lines 93-94: Check GPU upload path when hardware is available; implemented by `if gpuDeviceCount>0`.

### Control flow inferred from the code

- Line 94: conditional branch on `gpuDeviceCount>0`.

### Key state/data transformations

- Lines 8: computes `a` using `a=randn(2,2)+1i*randn(2,2)`.
- Lines 9: computes `b` using `b=randn(3,3)+1i*randn(3,3)`.
- Lines 10: computes `c` using `c=sprandn(2,2,0.75)+1i*sprandn(2,2,0.75)`.
- Lines 11: computes `d` using `d=randn(3,3)+1i*randn(3,3)`.
- Lines 12: computes `left` using `left=randn(5,6)+1i*randn(5,6)`.
- Lines 13: computes `right` using `right=randn(6,4)+1i*randn(6,4)`.
- Lines 14: computes `pre` using `pre=randn(8,6)+1i*randn(8,6)`.
- Lines 15: computes `suf` using `suf=randn(6,7)+1i*randn(6,7)`.
- Lines 18: computes `p` using `p=polyadic({{a,b},{c,d}})`.
- Lines 19: computes `ref` using `ref=kron(a,b)+kron(c,d)`.
- Lines 28: computes `p_pref` using `p_pref=prefix(left,p)`.
- Lines 30: computes `ref_pref` using `ref_pref=left*ref*right`.
- Lines 37: computes `q` using `q=polyadic({{d.',a.'}})`.
- Lines 38: computes `ref_q` using `ref_q=kron(d.',a.')`.
- Lines 51: computes `r` using `r=polyadic({{randn(2,4)+1i*randn(2,4),randn(3,2)+1i*randn(3,2)}})`.
- Lines 52: computes `ref_r` using `ref_r=kron(r.cores{1}{1},r.cores{1}{2})`.
- Lines 57: computes `k_num` using `k_num=randn(2,2)+1i*randn(2,2)`.
- Lines 70: computes `p_nan` using `p_nan=prefix(NaN,p)`.

### Local helper functions

- Line 104: `assert_small()` — `function assert_small(value,threshold)`.
  - Representative operation: `if value>=threshold`.
  - Representative operation: `error('numeric tolerance test failed.')`.

## Implementation structure

- Unit tests for advanced polyadic functionality.
- Get random test matrices
- Build a reference polyadic and its matrix form
- Check constructor, full, inflate, and validate
- Check prefixes, suffixes, size, and emptiness
- Check addition and subtraction paths
- Check multiplication paths
- Check Kronecker products
- Check transpose operations
- Check finiteness and internal non-zero counts
- Check zero-dimension behaviour
- Check nested simplification paths

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `sprandn()`, `polyadic()`, `assert_small()`, `inflate()`, `validate()`, `prefix()`, `suffix()`, `assert()`, `isequal()`, `transpose()`, `ctranspose()`, `allfinite()`, `all()`, `ref()`, `nnz()`, `simplify()`.
