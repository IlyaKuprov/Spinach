# kernel/optimcon/cubic_interp.m

Source: [kernel/optimcon/cubic_interp.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/optimcon/cubic_interp.m)
Wiki: [Spinach documentation for cubic_interp.m](https://spindynamics.org/wiki/index.php?title=cubic_interp.m)

## Purpose

Maximises the cubic Hermite interpolant defined by two function values and their directional derivatives. The interpolation interval has boundaries <code>end_a</code> and <code>end_b</code>; the two anchor points are <code>alpha_a</code> and <code>alpha_b</code>, with values <code>f_a</code> and <code>f_b</code> and directional derivatives <code>dir_deriv_a</code> and <code>dir_deriv_b</code>.

The routine normalises the anchor coordinate using <code>s = (alpha-alpha_a)/(alpha_b-alpha_a)</code>, constructs the cubic from the two endpoint values and derivatives, and evaluates its real stationary points that lie within the normalised interval together with both interval boundaries. It returns the point with the largest interpolated value, transformed back to the <code>alpha</code> coordinate. The second output is the value of the cubic model at that point, not a new evaluation of the underlying objective. This helper returns no gradient or adjoint.

## Call and validated inputs

<code>[alpha,fx] = cubic_interp(end_a,end_b,alpha_a,alpha_b,f_a,dir_deriv_a,f_b,dir_deriv_b)</code>

All eight inputs must be finite, real, numeric scalars, and <code>alpha_a</code> must differ from <code>alpha_b</code>. The code does not require <code>end_a</code> and <code>end_b</code> to differ, nor does it require the anchor points to lie inside the interpolation interval. There are no default input values.
