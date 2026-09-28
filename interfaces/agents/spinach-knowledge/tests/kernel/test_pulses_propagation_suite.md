# tests/kernel/test_pulses_propagation_suite.m

- Signature: `result=test_pulses_propagation_suite()`
- Output: regression result with explanatory messages.

## Checks
- RF polar–Cartesian gradient/Hessian round-trips for `f=sum(r.^2)=sum(x.^2+y.^2)`: Cartesian gradient `(2*x,2*y)`, diagonal Hessian blocks `2*I`, mixed blocks zero.
- `isergen`: second order `(HL+HR)/2+(1i*dt/6)*[HL,HR]`; fourth order `(HL+4*HM+HR)/6+(1i*dt/12)*[HL,HR]`.
- Constant-generator `PWCL`, `LG2`, `LG4`, `RKMK4`, and `LG4A` match `step`.
- `rsequence(1,4,1,1,1000,...)` phases `[base;-base]`, `base=[pi/4;-pi/4;pi/4;-pi/4]`; amplitude `4000*pi`, duration `1/4000`. Zero-RF compiler index map `[1;2;1]` gives identity propagators.
