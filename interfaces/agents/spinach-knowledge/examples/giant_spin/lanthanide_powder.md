# examples/giant_spin/lanthanide_powder.m

- Signature: `lanthanide_powder()`

## Purpose

Powder spectrum of Gd(III) with ZFS up to 4th spherical rank using the giant spin Hamiltonian formalism in a sweepable 400 MHz NMR magnet and microwaves at 263.2 GHz. Odd ranks are zero because a zero-field Hamiltonian must be even under time reversal. The 4th rank terms are converted from the Stevens parameters b40=4e-4 cm^-1 and b44=-2e-4 cm^-1 reported for Gd(III) in tetragonal BaTiO3 by Rimai and deMars (https://doi.org/10.1103/PhysRev.127.702). Calculation time: seconds.

## Physical / mathematical content

This is a single effective Gd(III) giant spin (E8) with isotropic Zeeman scalar 1.9918 and giant-spin terms through rank 4. The source sets the rank-1 and rank-3 coefficient arrays to zero; the rank-2 array is [0, 0, -4.65e8, 0, 0], and the rank-4 array is [-2.00e5, 0, 0, 0, 3.34e6, 0, 0, 0, -2.00e5]. All four giant-spin Euler-angle arrays are [0, 0, 0].

## Numerical / algorithmic content

The calculation uses a 1 T magnet and the zeeman-hilb basis with no approximation. It samples orientations on rep_2ang_100pts_sph, uses 263.2 GHz microwaves, linewidth parameter 2e-4, int_tol 10.0 and tm_tol 0.1, and sweeps 4096 points over [9.32, 9.56] T with rspt_order set to Inf. The initial state is the negative E8 Lz operator; fieldsweep computes the spectrum, which is plotted against its returned field axis.
