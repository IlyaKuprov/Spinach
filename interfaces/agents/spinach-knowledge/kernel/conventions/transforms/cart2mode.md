# kernel/conventions/transforms/cart2mode.m

[Source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/conventions/transforms/cart2mode.m) · [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=cart2mode.m)

## Contract

cart2mode projects Cartesian first- or second-derivative data for an interaction parameter onto one or two mass-weighted normal modes, including the mode's zero-point displacement scale. It produces derivatives with respect to the dimensionless coordinate (a+a')/sqrt(2) used by the bosonic mode interface in create.m. It is a derivative-coordinate transform, not an orientation-grid or spatial-coordinate transform.

## Inputs and dimensions

Let N be the number of atoms. Cartesian degrees of freedom are ordered [x1 y1 z1 x2 y2 z2 ...].

- cart_derivs: first derivatives in Hz per Angstrom, shape [d1 d2 3N]; or second derivatives in Hz per Angstrom squared, shape [d1 d2 3N 3N].
- eigvecs: documented as orthonormal mass-weighted normal-mode eigenvectors. The first-order case uses one [3N 1] vector; the second-order case uses two columns [3N 2].
- masses: positive atomic masses in unified atomic mass units, as an [N 1] column.
- frqs: positive mode frequency in Hz, a scalar for first order or a [1 2] vector for second order.

Each mode's Cartesian displacement scale is based on the zero-point amplitude sqrt(hbar/(m*omega)), with omega=2*pi*frqs; the implementation converts this displacement to Angstrom and converts masses from unified atomic mass units to kg. It contracts each first derivative once with its mode scale, or each second derivative with both scales.

## Output and use

mode_derivs is a [d1 d2] array in Hz. It is intended for the corresponding inter.modes.coupling_mod cell (d1=3, d2=3) or inter.modes.zeeman_mod cell (d1=1, d2=3). These dimensions describe the interaction parameter data, not atoms or modes. The routine returns raw Taylor derivatives; Spinach applies the Taylor-series one-half factor internally for second-order terms. Derivative data in wavenumbers or meV must first be converted to Hz with icm2hz.m or mev2hz.m. Zero and negative frequencies are rejected because the zero-point scaling is undefined.

## Source-supported use

The documented call is mode_derivs=cart2mode(cart_derivs,eigvecs,masses,frqs). The input descriptions above are the examples of supported first- and second-order array layouts given by the source; no concrete numerical dataset is supplied.
