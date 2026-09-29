# etc/textbook/rlx_nqi.m

## Use

[r1,r2,t1,t2]=rlx_nqi(I,omega,C_q,eta_q,tau_c) evaluates Redfield quadrupolar relaxation for a nucleus in an isotropically tumbling liquid.

## Inputs

- I: nuclear spin quantum number; the implementation requires integer or half-integer I >= 1.
- omega: nuclear Zeeman angular frequency in rad/s.
- C_q: quadrupolar coupling constant e^2*q*Q/h in Hz.
- eta_q: quadrupolar tensor asymmetry parameter. The source does not specify a unit or impose a further range.
- tau_c: positive rotational correlation time in seconds.

All five inputs must be real numeric scalars. The code does not state a finiteness check for these scalar inputs.

## Calculation and outputs

The routine constructs the quadrupolar tensor with eeqq2nqi, converts its Hz coupling to angular units with 2*pi, and uses the tensor's rank-2 Blicharski invariant. With isotropic rotational diffusion D = 1/(6*tau_c), r1 samples rank-2 spectral densities at omega and 2*omega; r2 uses zero, omega, and 2*omega. It returns reciprocal relaxation rates as times:

- r1: longitudinal relaxation rate, Hz; t1 = 1/r1: longitudinal relaxation time, seconds.
- r2: transverse relaxation rate, Hz; t2 = 1/r2: transverse relaxation time, seconds.

## Scope and source

This is the isotropic-tumbling quadrupolar model; omega is supplied directly rather than derived from a field. The source gives no allowed interval for eta_q, so none is asserted here. Source: [implementation](https://github.com/IlyaKuprov/Spinach/blob/main/etc/textbook/rlx_nqi.m) · [Spin Dynamics Wiki](https://spindynamics.org/wiki/index.php?title=rlx_nqi.m).
