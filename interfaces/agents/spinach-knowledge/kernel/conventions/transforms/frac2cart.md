# kernel/conventions/transforms/frac2cart.m

- Signature: `[XYZ,va,vb,vc]=frac2cart(a,b,c,alp,bet,gam,ABC)`

## Purpose

Converts fractional crystallographic coordinates to Cartesian coordinates for a unit cell described by three lengths and three angles.

## Parameters / inputs

- `a,b,c`: positive real scalar unit-cell dimensions.
- `alp,bet,gam`: real scalar unit-cell angles in degrees.
- `ABC`: real `N x 3` array of fractional coordinates.

## Outputs

- `XYZ`: `N x 3` Cartesian coordinates, in units consistent with the unit-cell dimensions.
- `va,vb,vc`: primitive lattice vectors, returned as the columns of the cell transformation matrix.

## Method

The function constructs the triclinic-cell transformation matrix from the lengths and angle cosines/sines, applies it to each coordinate row, and returns its three columns as the lattice vectors.

Source: [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=frac2cart.m)
