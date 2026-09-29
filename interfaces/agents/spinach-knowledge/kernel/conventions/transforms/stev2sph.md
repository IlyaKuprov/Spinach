# kernel/conventions/transforms/stev2sph.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/conventions/transforms/stev2sph.m) · [Spin Dynamics Wiki page](https://spindynamics.org/wiki/index.php?title=stev2sph.m)

## Purpose and convention

Converts coefficients of Stevens operators into coefficients of irreducible spherical tensor operators. For each rank `k`, the input has `2*k+1` real coefficients ordered by projections `-k:k` (increasing); the output has `2*k+1` complex coefficients ordered `k:-1:-k` (decreasing). The scaling factors are dimensionless: this conversion does not change the coefficients’ physical units.

## Conversion

For ranks `k=1,...,12`, the implementation builds a diagonal-plus-antidiagonal matrix `A=diag(u)+fliplr(diag(v))` and returns `transpose(B_in'*A)`, where `B_in` is the input column. Its diagonal and antidiagonal entries are:

- `u=[-1i*(-1).^(k:-1:1)'; 1; ones(k,1)].*a{k}`
- `v=[1i*ones(k,1); 0; (-1).^(1:k)'].*a{k}`

The rank-specific scale vector `a{k}` is defined in the source. Its squared entry for projection magnitude `q>0` is `2^(k-2)*P(k,q)/C(k,q)^2`; for `q=0` it is `2^k*P(k,0)/C(k,0)^2`. Here `P(k,q)=product((k+p)*(k-p+1), p=q+1,...,k)`, and `C(k,q)` is the corresponding Stevens coefficient from [`stevens.m`](../../operators/stevens.md), doubled when `k` is even and `q` is odd. The source's explicit rank tables implement these factors through rank 12. The formula reproduces the published ranks 1–6; the cited published source is [J. Phys. C: Solid State Phys. 18, 1429 (1985)](https://doi.org/10.1088/0022-3719/18/7/009). The rank-7–12 rational factors are based on the integer coefficient table of Ryabov, *J. Magn. Reson.* 140, 141 (1999); the source notes large-prime denominators at rank 9, projections 1 and 2.

The input rank must be a finite, real numeric scalar integer from 1 through 12. `B_in` must be a finite, real numeric column with exactly `2*k+1` elements. The result is a complex column vector of the same length, suitable as coefficients for [`irr_sph_ten.m`](../../operators/irr_sph_ten.md).
