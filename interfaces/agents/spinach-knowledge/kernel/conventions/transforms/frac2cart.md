# kernel/conventions/transforms/frac2cart.m

**MATLAB source:** [kernel/conventions/transforms/frac2cart.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/conventions/transforms/frac2cart.m)
**Spinach Wiki:** [frac2cart.m](https://spindynamics.org/wiki/index.php?title=frac2cart.m)

## Purpose and units

Converts rows of fractional crystallographic coordinates to Cartesian coordinate rows using the transformation matrix computed from the unit-cell lengths and angles. The source expects positive unit-cell lengths `a,b,c` and angles `alp,bet,gam` in degrees. Its guard checks that each length is a numeric real scalar and rejects values `<= 0`, but does not test finiteness (so a `NaN` value is not screened by that comparison). It does not specify a particular length unit: if the three lengths use a chosen unit, Cartesian coordinates and lattice vectors have that same unit.

## Inputs and outputs

Signature: `[XYZ,va,vb,vc]=frac2cart(a,b,c,alp,bet,gam,ABC)`

- `a,b,c`: numeric real scalars, each rejected if `<= 0`; no separate finite-value check is made.
- `alp,bet,gam`: numeric real scalars in degrees.
- `ABC`: numeric real `N x 3` array of fractional coordinates.
- `XYZ`: `N x 3` Cartesian coordinate array.
- `va,vb,vc`: three-column lattice vectors (each is `3 x 1`).

Validation does not require finite angles or coordinates and does not impose an angle range or a separate cell-geometry check; the values are used directly by the trigonometric formulas below.

## Transformation

The source computes the scalar

`v = a*b*c*sqrt(1-cosd(alp)^2-cosd(bet)^2-cosd(gam)^2+2*cosd(alp)*cosd(bet)*cosd(gam))`

and the matrix

`T = [a, b*cosd(gam), c*cosd(bet); 0, b*sind(gam), c*(cosd(alp)-cosd(bet)*cosd(gam))/sind(gam); 0, 0, v/(a*b*sind(gam))]`.

For each input row, the implemented mapping is `XYZ = (T*ABC')'`, equivalently `ABC*T'`. The primitive vectors are the columns of `T`: `va=T(:,1)`, `vb=T(:,2)`, and `vc=T(:,3)`.
