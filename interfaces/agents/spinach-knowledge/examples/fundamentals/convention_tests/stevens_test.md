# examples/fundamentals/convention_tests/stevens_test.m

- MATLAB implementation: [examples/fundamentals/convention_tests/stevens_test.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/fundamentals/convention_tests/stevens_test.m)

- Signature: `stevens_test()`

## Question tested

For spin quantum number 12 (multiplicity 25), the test compares Spinach's Stevens matrices with explicit spin-operator expressions for every component `q=-k,…,k` at ranks `k=6,4,2`. It also checks that the five rank-2 Stevens matrices match the rank-2 irreducible spherical tensor (IST) combinations used in the source. The comment describes the explicit expressions as being from the literature, but supplies no bibliographic citation.

## Explicit-operator comparison

The source obtains the Spinach matrices from `stevens(mult,k,q)` and builds the comparison matrices from `S=pauli(mult)` with `s=sqn*(sqn+1)`, using polynomials in `S.z`, powers of `S.p` and `S.m`, and symmetrised products through `acomm`. For example, at rank 2 it uses `O_2,0 = 3*S.z^2 - s*S.u` and `O_2,+2 = (S.p^2+S.m^2)/2`; the negative-q expression has the source's imaginary phase, `O_2,-2 = (S.p^2-S.m^2)/(2i)`. Thus the comparisons exercise both normalisation and signed-component/phase conventions, rather than only the diagonal component.

For each rank the code forms a vector of individual matrix 1-norm residuals, then applies MATLAB's default vector norm to that vector. The error thresholds are:

- rank 6: `norm(diffs) > 1e-4`;
- rank 4: `norm(diffs) > 1e-7`;
- rank 2: `norm(diffs) > 1e-12`.

Here `norm(diffs)` is the default 2-norm of the residual vector; it is not the vector's 1-norm.

## Rank-2 IST comparison

With `T=irr_sph_ten(mult,2)`, the source defines `T2p2=T{1}`, `T2p1=T{2}`, `T20=T{3}`, `T2m1=T{4}`, and `T2m2=T{5}`. Its matching combinations are `O_2,+2=T2m2+T2p2`, `O_2,+1=(T2m1-T2p1)/2`, `O_2,0=sqrt(6)*T20`, `O_2,-1=i*(T2m1+T2p1)/2`, and `O_2,-2=i*(T2m2-T2p2)`. The norm of the vector of five matrix 1-norm residuals must not exceed `1e-10` or the source errors.

The page records the implemented comparisons and thresholds; it does not assert that the test was executed or passed.
