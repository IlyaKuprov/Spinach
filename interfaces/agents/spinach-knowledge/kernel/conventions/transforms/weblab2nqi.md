# kernel/conventions/transforms/weblab2nqi.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/conventions/transforms/weblab2nqi.m) · [Spin Dynamics Wiki page](https://spindynamics.org/wiki/index.php?title=weblab2nqi.m)

## Purpose and inputs

Converts the Weblab one-cone model parameters ([diagram](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/conventions/weblab_cone.png)) into quadrupolar coupling tensors for multiple sites. `C_q` is `e^2*q*Q/h` in Hz; `eta_q` is the dimensionless quadrupolar asymmetry; `I` is the spin quantum number; and `alpha`, `theta`, and, where accepted, `phi` are model angles in radians. The generated site tensors are 3x3 matrices in Hz.

The inputs must be numeric, real scalars. `I` must also satisfy `I>=1` and `2*I+1` integer (integer or half-integer spins starting at 1). The code does not impose finite-value checks or an interval on `eta_q` or the angles.

## Output modes and azimuths

The number and arrangement of outputs determines the accepted input form:

- Two outputs require six inputs including `phi`; the site azimuths are `-phi/2` and `+phi/2`.
- Three outputs require six inputs; the azimuths are `-phi`, `0`, and `+phi`.
- Four outputs require five inputs, with `phi` omitted; the fixed azimuths are `0`, `pi/2`, `pi`, and `3*pi/2`.
- Six outputs require five inputs, with `phi` omitted; the fixed azimuths are `0`, `pi/3`, `2*pi/3`, `pi`, `4*pi/3`, and `5*pi/3`.

For every site, the wrapper calls `eeqq2nqi(C_q,eta_q,I,[azimuth theta alpha])`. That routine uses the ZYZ active Euler convention in radians. Its principal values are `XX=-C_q*(1-eta_q)/(4*I*(2*I-1))`, `YY=-C_q*(1+eta_q)/(4*I*(2*I-1))`, and `ZZ=C_q/(2*I*(2*I-1))`; it rotates `diag([XX YY ZZ])` by `R=euler2dcm([azimuth theta alpha])` as `Q=R*diag([XX YY ZZ])*R'`, then removes the trace-rounding component and symmetrises the matrix. Each output `Qn` is one such 3x3 Hz tensor.

Only 2, 3, 4, or 6 outputs are supported. Two- and three-output calls require exactly six inputs; four- and six-output calls require exactly five and reject a supplied `phi`. `phi` must be scalar in the six-input modes. The five-input modes set it internally to an empty array.
