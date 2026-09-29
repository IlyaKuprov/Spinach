# tests/kernel/test_lowlevel_utilities_suite.m

## Purpose

Regression test suite for cheap deterministic low-level utility functions in the Spinach kernel. It checks small numerical helpers, matrix filters, integer type selection, and analytic line-shape definitions against explicit answers.

## Behaviour

The function announces the test target with `fprintf('TESTING: Low-level utility functions\n')` and initialises a test result object via `new_test_result('kernel/lowlevel_utilities_suite', 'Low-level utility functions', 'cheap scalar and matrix helpers must preserve their documented algebraic definitions.')`.

Using test matrices `A=[1 2;3 4]`, `B=[0 1;-1 2]`, `C=[2 1;0 4]`, `D=[1+1i 2-1i;3 4i]`, `E=[2 0;1i -3]`, and `M=reshape(1:16,4,4)`, the suite verifies, each with tolerances `1e-15` (or `1e-14` for `keep_rank`):

- `comm(A,B)` equals `A*B-B*A`.
- `rocomm({A,B,A})` equals `comm(comm(A,B),A)`, nesting commutators left to right over the supplied cell array.
- `remtrace(C)` equals `[-1 1;0 1]`, subtracting `trace(A)/dim` times the unit matrix.
- `remncomm(H,eye(2),[1;2])` with `H=[2 1+2i;1-2i 5]` returns `diag(diag(H))`: in the eigenbasis of a non-degenerate `B` only the diagonal part commutes with `B`.
- `remncomm(H,eye(3),[1;1;2])` with `H=[2 1+2i 3;1-2i 5 1i;3 -1i 1]` returns `[2 1+2i 0;1-2i 5 0;0 0 1]`: in the eigenbasis of a degenerate `B` the whole degenerate block of `A` commutes with `B`.
- `remncomm(H,eye(2),[1e12;1e12+1])` returns `diag(diag(H))`: eigenvalues of `B` differing by one on a large common offset are not degenerate.
- `remncomm(H,eye(3),[0;1;1e12])` returns `diag(diag(H))`: a small splitting is not merged by a wide spectral span.
- `remncomm(H,eye(2),[1e12;1e12+0.01])` returns `diag(diag(H))`: the degeneracy criterion does not depend on the energy origin.
- `hdot(D,E)` equals `trace(D'*E)`, the Frobenius inner product.
- `atranspose(A)` equals `[4 2;3 1]`, reflecting a matrix across the anti-diagonal.
- `killcross(M,[2 4],[1 3])` equals `[0 0 0 0;2 0 10 0;0 0 0 0;4 0 12 0]`, zeroing the requested columns and rows.
- `killdiag(M,1)` equals `[0 5 9 13;2 0 10 14;3 7 0 15;4 8 12 0]`, zeroing the requested diagonal band.
- `keep_rank(S,2)` with `S=diag([5 2 1])` equals `diag([5 2 0])`, reconstructing from the requested leading singular values (tolerance `1e-14`).
- `frob_chop([5 4 0.3 0.1],0.5)` equals `2`, keeping the smallest rank whose discarded tail is below tolerance.

Analytic line shapes are checked at `x=[0 1]` with `g_fwhm=2*sqrt(2*log(2))`:

- `gaussfun(x,g_fwhm)` equals `exp(-(x.^2)/2)/sqrt(2*pi)`: with sigma one, `gaussfun` is the standard normal density.
- `[lor_r,lor_i]=lorentzfun(0,2*pi,2,x,0)` gives `lor_r=[1 1/2]` (zero-phase Lorentzian real part is `1/(1+x^2)` for the chosen parameters) and `lor_i=[0 1/2]` (imaginary part is `x/(1+x^2)`).
- `spden(2,10,0)` equals `1/300`: rank-two rotational spectral density at zero frequency is `tau_c/(2L+1)`.
- `spden(2,10,60)` equals `1/600`: when `omega*tau_c` is one, the Lorentzian spectral density is halved.

Integer type selection is checked at promotion boundaries using `test_true`:

- `min_int_type(127,'signed')` returns `'int8'` (127 is representable in signed 8-bit storage).
- `min_int_type(128,'signed')` returns `'int16'` (128 requires signed 16-bit storage).
- `min_int_type(255,'unsigned')` returns `'uint8'` (255 is representable in unsigned 8-bit storage).
- `min_int_type(256,'unsigned')` returns `'uint16'` (256 requires unsigned 16-bit storage).

## Inputs and outputs

- **result** (output) - regression test result with explanatory messages, returned by the function.

The function takes no inputs.

## References

- Source: [tests/kernel/test_lowlevel_utilities_suite.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_lowlevel_utilities_suite.m)
