# kernel/conventions/transforms/cart2mode.m

- Signature: `mode_derivs=cart2mode(cart_derivs,eigvecs,masses,frqs)`

## Purpose

Converts Cartesian derivatives of spin Hamiltonian parameters into derivatives with respect to dimensionless normal-mode coordinates for `inter.modes.coupling_mod` and `inter.modes.zeeman_mod` in the bosonic mode specification used by `create.m`.
## Physical / mathematical content

The dimensionless coordinate of a mode is `(a+a')/sqrt(2)`. For a Cartesian degree of freedom with mass `m`, its displacement per unit mode coordinate is the mass-weighted eigenvector component times `sqrt(hbar/(m*omega))`, where `omega=2*pi*frqs`. First derivatives are contracted with one displacement vector; second derivatives are contracted with two, one for each mode.
## Numerical / algorithmic content

Atomic masses are expanded over their three Cartesian degrees of freedom and converted from unified atomic mass units to kilograms. The code computes displacement scales in Angstrom, then sums the Cartesian derivative array against one scale vector for first order or two scale vectors for second order. It validates finite real inputs, positive masses and frequencies, matching dimensions, and unit-normalised eigenvector columns.
## Parameters / inputs

- cart_derivs -first derivatives of an interaction with
- respect to Cartesian displacements, in Hz
- per Angstrom, an array of dimension
- [d1 d2 3N] where N is the number of atoms
- and the Cartesian degrees of freedom are
- ordered [x1 y1 z1 x2 y2 z2 ...]; or second
- derivatives in Hz per Angstrom squared, an
- array of dimension [d1 d2 3N 3N]
- eigvecs -orthonormal mass-weighted normal mode
- eigenvector, a [3N 1] column vector for
- the first order case; two such vectors
- as a [3N 2] array for the second order
- case
- masses -atomic masses in unified atomic mass
- units, an [N 1] column vector
- frqs -mode frequency in Hz, a positive scalar
- for the first order case; a [1 2] vector
- for the second order case

## Outputs

- mode_derivs -derivatives with respect to dimensionless
- mode coordinates, in Hz, a [d1 d2] array
- to be placed into the corresponding cell
- of inter.modes.coupling_mod (d1=3, d2=3)
- or inter.modes.zeeman_mod (d1=1, d2=3)
- Note: raw Taylor derivatives are returned; the 1/2 factors of
- the Taylor expansion are applied by Spinach internally.
- Derivative data in wavenumbers or meV should be conver-
- ted into Hz with icm2hz.m or mev2hz.m beforehand. Zero
- and negative frequency modes are rejected because their
- zero-point scaling is undefined.

## Implementation structure

The main function calls the local `grumble` validator, constructs the zero-point displacement scale vectors, and contracts them with `cart_derivs`. The local validator enforces the distinct eigenvector and frequency shapes required for first- and second-order conversion.
