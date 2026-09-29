# examples/fundamentals/convention_tests/spsk_test.m

- MATLAB implementation: [examples/fundamentals/convention_tests/spsk_test.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/fundamentals/convention_tests/spsk_test.m)

- Signature: `spsk_test()`

## Question tested

The test checks the span–skew convention used by `spsk2mat`: do the isotropic value, span, skew, and Euler angles reconstruct the same rotated tensor as direct rotation of its principal-value diagonal matrix?

## Construction and criterion

It draws `xx=rand()`, `yy=rand()+3`, and `zz=rand()+6`, so the generated values are ordered and distinct, then draws Euler angles `alp` and `gam` from `[0,2π)` and `bet` from `[0,π)`. With `R=euler2dcm(alp,bet,gam)`, the direct matrix is `AM=R*diag([xx yy zz])*R'`. The convention parameters in the source are

- `iso=(xx+yy+zz)/3`;
- `sp=zz-xx`;
- `sk=3*(yy-iso)/sp`.

The Spinach construction is `AS=spsk2mat(iso,sp,sk,alp,bet,gam)`. It tests whether `norm(AM-AS,1)<1e-6`; the source displays its success message only on that branch and otherwise raises an error. This checks the parameterisation for generated diagonal eigenvalues and orientations, not every possible tensor input.
