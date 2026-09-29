# kernel/utilities/spher_harmon.m

## Purpose

Evaluates spherical harmonics Y(l,m,theta,phi) at user-specified polar and azimuthal angles.

## Behaviour

- Syntax: `Y=spher_harmon(l,m,theta,phi)`.
- Input consistency is enforced by an internal `grumble` subfunction, which errors out when: `l` is not a non-negative real integer; `m` is not a real integer in the interval `[-l,l]`; `theta` or `phi` is not numeric and real.
- Schmidt-normalised associated Legendre functions are obtained with MATLAB's `legendre(l,cos(theta),'sch')`, reshaped to `[l+1 numel(theta)]`, and the row `abs(m)+1` is extracted and reshaped back to the size of `theta`.
- The spherical harmonic is assembled as `sqrt((2*l+1)/(4*pi))*S.*exp(1i*m*phi)`; for nonzero `m` an additional division by `sqrt(2)` is applied.
- If `m>0` and `m` is odd, the sign of `Y` is flipped.

## Inputs and outputs

Inputs:

- `l` — L quantum number; non-negative real integer.
- `m` — M quantum number; real integer from `[-l,l]`.
- `theta` — array of theta angles in radians; numeric and real.
- `phi` — array of phi angles in radians; numeric and real.

Output:

- `Y` — array of spherical harmonics evaluated at the specified angles.

## References

- Source: <https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/spher_harmon.m>
- Spin Dynamics Wiki: <https://spindynamics.org/wiki/index.php?title=spher_harmon.m>
