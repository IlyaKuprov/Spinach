# kernel/conventions/transforms/qform2sph.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/conventions/transforms/qform2sph.m) · [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=qform2sph.m)

- Signature: `[r0,r1,r2]=qform2sph(A)`

## Behaviour

For the normalised quadratic form `[x y z]*A*[x y z]'/norm([x y z],2)^2`, the returned coefficients satisfy `quadratic_form = sum(r_LM*Y_LM)`. The input is a real symmetric 3x3 matrix.

- `r0 = (2/3)*sqrt(pi)*trace(A)`.
- `r1 = [0 0 0]`.
- `r2` has five entries ordered by decreasing m: m=2, 1, 0, -1, -2. With `c=sqrt(2*pi/15)`, they are:
  - `r2(1)=c*(A(1,1)-A(2,2)-1i*(A(1,2)+A(2,1)))`
  - `r2(2)=-c*(A(1,3)+A(3,1)-1i*(A(2,3)+A(3,2)))`
  - `r2(3)=-(2/3)*sqrt(pi/5)*(-2*A(3,3)+A(2,2)+A(1,1))`
  - `r2(4)=c*(A(1,3)+A(3,1)+1i*(A(2,3)+A(3,2)))`
  - `r2(5)=c*(A(1,1)-A(2,2)+1i*(A(1,2)+A(2,1)))`

No unit conversion is applied; the coefficient scale follows the units of `A`.

## Input and outputs

The explicit input check requires `A` to be numeric, real, symmetric, and exactly 3x3. Outputs are scalar `r0`, a 1x3 zero row vector `r1`, and a 1x5 row vector `r2`.
