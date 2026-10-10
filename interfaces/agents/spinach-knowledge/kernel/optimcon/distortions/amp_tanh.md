# kernel/optimcon/distortions/amp_tanh.m

- Signature: [w,J]=amp_tanh(w,sat_lvls)
- MATLAB source: [amp_tanh.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/optimcon/distortions/amp_tanh.m)

## Purpose

Models radial amplifier compression independently for each in-phase/quadrature (XY) channel pair and time slice. For a pair with components X and Y, let r=sqrt(X^2+Y^2) and let a be its saturation level. The output has radial amplitude a*tanh(r/a) and retains the input direction; at r=0 the pair is unchanged.

## Inputs and units

- w is a real numeric waveform in rad/s nutation-frequency units. Columns are time slices; rows are ordered XYXY..., one X/Y pair per control channel. The number of rows must be even.
- sat_lvls is a numeric, real array with one finite, strictly positive value per XY pair. Each value is the limiting radial output amplitude, in the same amplitude units as w.

## Outputs and derivative

- w is the distorted waveform, with the input shape and units.
- J is optional. When requested, it is a sparse Jacobian of the vectorised output with respect to the vectorised input. For r>0, set s=a*tanh(r/a)/r and c=(1/cosh(r/a)^2-s)/r^2; the local XY block is J_pair=[s+c*X^2,c*X*Y;c*X*Y,s+c*Y^2]. At r=0 the implementation uses s=1 and c=0, giving the identity block. There are no cross-time or cross-channel terms. GPU-computed block values are gathered to host memory before sparse assembly.

## References

- Spinach documentation: https://spindynamics.org/wiki/index.php?title=amp_tanh.m
