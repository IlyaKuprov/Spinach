# etc/textbook/r2csa2tauc.m

- MATLAB implementation: [etc/textbook/r2csa2tauc.m](https://github.com/IlyaKuprov/Spinach/blob/main/etc/textbook/r2csa2tauc.m)

**Signature:** `tauc=r2csa2tauc(R2,del_sq,B0,isotope)`

## Purpose

Finds the three algebraic rotational-correlation-time solutions of the cubic relation implemented for a transverse CSA relaxation rate. The executable signature and parameter description use `R2`; the prototype in the source comment alone says `R1`.

## Inputs

- `R2` — positive real numeric scalar; the source documents the transverse relaxation rate in Hz.
- `del_sq` — positive real numeric scalar, documented as the second-rank CSA invariant (see `blinv.m`). It is dimensionless squared fractional shielding: `rlx_csa` obtains it as the second output of `blinv(1e-6*Z_ppm)`, not from the raw ppm tensor.
- `B0` — real numeric scalar magnetic field in tesla. The consistency check does not require it to be positive.
- `isotope` — character array identifying the isotope, for example `'1H'`.

## Calculation and output

The code obtains `omega=-B0*spin(isotope)` and evaluates three closed-form roots of the cubic (the source identifies the expressions as Mathematica-generated). The output `tauc` is a three-element vector of algebraic solutions in seconds; the function does not select a single positive/physical root. For the roots, it raises an error only when all three entries fail the real-valued check, with the message “no real solutions - physically impossible case.”

## Reference

[Spinach Wiki: r2csa2tauc.m](https://spindynamics.org/wiki/index.php?title=r2csa2tauc.m)
