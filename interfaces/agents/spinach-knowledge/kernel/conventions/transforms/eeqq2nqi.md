# kernel/conventions/transforms/eeqq2nqi.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/conventions/transforms/eeqq2nqi.m) · [Spin Dynamics Wiki](https://spindynamics.org/wiki/index.php?title=eeqq2nqi.m)

## Conversion

This function maps the quadrupolar specification C_q, eta_q, and I to a 3x3 coupling tensor Q in Hz. C_q = e^2*q*Q/h is the quadrupolar coupling constant in Hz; eta_q is the quadrupolar tensor asymmetry parameter; I is the spin quantum number. The source gives the principal-axis components:

~~~text
XX = -C_q*(1-eta_q)/(4*I*(2*I-1))
YY = -C_q*(1+eta_q)/(4*I*(2*I-1))
ZZ =  C_q/(2*I*(2*I-1))
~~~

The three Euler angles eulers are in radians and specify the principal-axis frame orientation relative to the lab frame. The implementation calls euler2dcm(eulers) and applies its active ZYZ convention:

~~~text
R = euler2dcm(eulers)
Q = R*diag([XX YY ZZ])*R'
Q = Q - eye(3)*trace(Q)/3
Q = (Q+Q')/2
~~~

Thus Q is returned as a symmetric, traceless 3x3 tensor in Hz. The source notes that the denominator contains the spin quantum number squared, so the tensor falls off sharply with nuclear spin for the same anisotropy parameter.

## Inputs and constraints

C_q, eta_q, I, and eulers must be numeric and real. C_q, eta_q, and I must each be scalar; eulers must have exactly three elements (the validator does not require a particular row/column shape). I must be an integer or half-integer at least 1. The source imposes no additional range or finiteness check on C_q, eta_q, or the angles.
