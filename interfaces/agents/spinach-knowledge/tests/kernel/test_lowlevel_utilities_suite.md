# tests/kernel/test_lowlevel_utilities_suite.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_lowlevel_utilities_suite.m`
- Signature: `result=test_lowlevel_utilities_suite()`
- Total lines: 93

## Purpose

Tests cheap deterministic low-level utility functions. Syntax: result=test_lowlevel_utilities_suite()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 16-17: Announce the test target; implemented by `fprintf('TESTING: Low-level utility functions\n')`.
- Lines 19-22: State the utility target of the test; implemented by `result=new_test_result('kernel/lowlevel_utilities_suite', 'Low-level utility functions', 'cheap scalar and matrix helpers must preserve their documented algebraic defini…`.
- Lines 24-25: Define small test matrices; implemented by `A=[1 2;3 4]`.
- Lines 28-30: Check commutator and right-ordered nested commutator; implemented by `result=test_close(result,'comm',comm(A,B),A*B-B*A,1e-15,1e-15, 'the commutator is AB-BA')`.
- Lines 34-35: Check trace removal and commuting part extraction; implemented by `C=[2 1;0 4]`.
- Lines 42-43: Check Frobenius inner product and anti-diagonal transpose; implemented by `D=[1+1i 2-1i;3 4i]`.
- Lines 50-51: Check matrix wiping helpers; implemented by `M=reshape(1:16,4,4)`.
- Lines 59-60: Check rank and SVD truncation helpers; implemented by `S=diag([5 2 1])`.
- Lines 66-67: Check analytic line shapes and spectral density at points with closed-form values; implemented by `x=[0 1]`.
- Lines 81-83: Check minimum integer type selection at promotion boundaries; implemented by `result=test_true(result,'min_int_type signed int8',strcmp(min_int_type(127,'signed'),'int8'), '127 is representable in signed 8-bit storage')`.

### Key state/data transformations

- Lines 20-22: computes `result` using `result=new_test_result('kernel/lowlevel_utilities_suite', 'Low-level utility functions', 'cheap scalar and matrix helpers must preserve their documented algebraic defini…`.
- Lines 25: computes `A` using `A=[1 2;3 4]`.
- Lines 26: computes `B` using `B=[0 1;-1 2]`.
- Lines 35: computes `C` using `C=[2 1;0 4]`.
- Lines 38: computes `H` using `H=[2 1+2i;1-2i 5]`.
- Lines 43: computes `D` using `D=[1+1i 2-1i;3 4i]`.
- Lines 44: computes `E` using `E=[2 0;1i -3]`.
- Lines 51: computes `M` using `M=reshape(1:16,4,4)`.
- Lines 60: computes `S` using `S=diag([5 2 1])`.
- Lines 67: computes `x` using `x=[0 1]`.
- Lines 68: computes `g_fwhm` using `g_fwhm=2*sqrt(2*log(2))`.
- Lines 71: computes `[lor_r,lor_i]` using `[lor_r,lor_i]=lorentzfun(0,2*pi,2,x,0)`.

## Outputs

- result -regression test result with explanatory messages
- The test checks small numerical helpers, matrix filters, integer type
- selection, and analytic line-shape definitions against explicit answers.

## Implementation structure

- Tests cheap deterministic low-level utility functions. Syntax:
- result=test_lowlevel_utilities_suite()
- result -regression test result with explanatory messages
- The test checks small numerical helpers, matrix filters, integer type
- selection, and analytic line-shape definitions against explicit answers.
- Announce the test target
- State the utility target of the test
- Define small test matrices
- Check commutator and right-ordered nested commutator
- Check trace removal and commuting part extraction
- Check Frobenius inner product and anti-diagonal transpose
- Check matrix wiping helpers

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `test_close()`, `comm()`, `rocomm()`, `remtrace()`, `remncomm()`, `hdot()`, `atranspose()`, `killcross()`, `killdiag()`, `keep_rank()`, `frob_chop()`, `gaussfun()`, `lorentzfun()`, `spden()`, `test_true()`.
