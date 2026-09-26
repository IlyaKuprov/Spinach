# kernel/conventions/transforms/frac2cart.m

- Signature: `[XYZ,va,vb,vc]=frac2cart(a,b,c,alp,bet,gam,ABC)`

## Purpose

Converts fractional crystallographic coordinates to Cartesian coordinates. Syntax: [XYZ,va,vb,vc]=frac2cart(a,b,c,alpha,beta,gamma,ABC)

## Physical / mathematical content

- Convention and tensor-transform utilities. They convert among tensor parameterisations, coordinate systems, and unit systems; the underlying mathematics is linear algebra on rank-2 tensors and rotation representations.

## Numerical / algorithmic content

## Parameters / inputs

- a,b,c -three unit cell dimensions
- alp,bet,gam -three unit cell angles, degrees
- ABC -fractional atomic coordinates as
- Nx3 array of numbers

## Outputs

- XYZ -Cartesian atomic coordinates as
- Nx3 array of numbers
- va, vb, vc -primitive lattice vectors

## Implementation structure

- Converts fractional crystallographic coordinates to Cartesian
- coordinates. Syntax:
- [XYZ,va,vb,vc]=frac2cart(a,b,c,alpha,beta,gamma,ABC)
- a,b,c -three unit cell dimensions
- alp,bet,gam -three unit cell angles, degrees
- ABC -fractional atomic coordinates as
- Nx3 array of numbers
- XYZ -Cartesian atomic coordinates as
- va, vb, vc -primitive lattice vectors
- Check consistency
- Compute the transformation matrix
- Apply the transformation matrix
