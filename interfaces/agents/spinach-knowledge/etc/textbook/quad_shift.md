# etc/textbook/quad_shift.m

- MATLAB implementation: [etc/textbook/quad_shift.m](https://github.com/IlyaKuprov/Spinach/blob/main/etc/textbook/quad_shift.m)

- Signature: delta = quad_shift(Cq,eta,v0,S,m)

## Purpose

Calculates the second-order quadrupolar shift of the powder-pattern centre of gravity for the NMR transition |S,m> to |S,m-1>.

## Inputs

All five arguments are required; there are no defaults.

- Cq — real numeric scalar quadrupolar constant, in Hz.
- eta — real numeric scalar quadrupolar asymmetry parameter. The code does not enforce a range.
- v0 — real numeric scalar Larmor frequency, in Hz. The expression divides by v0, so use a nonzero frequency; the input check itself does not reject zero.
- S — real numeric scalar integer or half-integer spin, strictly greater than 1/2.
- m — real numeric scalar, integer or half-integer, satisfying 1-S <= m <= S. Although the error text describes an existing transition, the check does not require m to have the same integer/half-integer parity as S; ensure it is an actual spin-S projection for the intended physical transition.

## Expression and output

The implementation uses Samoson's Equation 3:

~~~matlab
delta=-1e6*(3/40)*(Cq/v0)^2*(1+eta^2/3)*...
           (S*(S+1)-9*m*(m-1)-3)/(S^2*(2*S-1)^2);
~~~

delta is the second-order shift in ppm. The source comment says this expression was checked against pure numerics; it does not specify the test setup.

## Source

Samoson, Equation 3, [Chemical Physics Letters (1985), DOI: 10.1016/0009-2614(85)85414-2](https://doi.org/10.1016/0009-2614(85)85414-2). [Spinach Wiki: quad_shift.m](https://spindynamics.org/wiki/index.php?title=quad_shift.m).
