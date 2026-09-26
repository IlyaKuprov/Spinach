# experiments/pseudocon/kpcs.m

- Signature: `[pcs_vals,pcs_cube]=kpcs(probden,chi,ranges,nxyz,method)`

## Purpose

Computes the three-dimensional distribution of pseudocontact shift field by solving Kuprov equation for PCS. Syntax: [pcs_cube,pcs_vals]=kpcs(probden,chi,ranges,nxyz)

## Physical / mathematical content

- Paramagnetic-pseudocontact inference routines. The mathematics includes inverse problems, tensor parameterisation, interpolation, and regularisation.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.
- Finite-difference discretisation appears in the implementation, so numerical accuracy depends on stencil order, boundary handling, and the balance between resolution and conditioning.

## Parameters / inputs

- probden -unpaired electron probability density (not spin density) cube
- chi -electron magnetic susceptibility tensor in cubic Angstroms
- ranges -Cartesian axis extents for the unpaired electron probability
- density cube as [xmin xmax ymin ymax zmin zmax] in Angstroms
- nxyz -nuclear coordinates as [x y z] with multiple rows) at which
- PCS is to be evaluated, in Angstroms
- Output:
- pcs_vals -pseudocontact shift in ppm at each nucleus
- pcs_cube -pseudocontact shift field on the same grid as the unpaired
- electron probability density supplied
- Note: minimal three-point schemes are used for the finite difference
- operators. Increase stencil size if you have enough memory.
- Note: for further information on the equations and algorithms used in this
- function see http://dx.doi.org/10.1039/C4CP03106G

## Implementation structure

- Computes the three-dimensional distribution of pseudocontact shift field
- by solving Kuprov equation for PCS. Syntax:
- [pcs_cube,pcs_vals]=kpcs(probden,chi,ranges,nxyz)
- probden -unpaired electron probability density (not spin density) cube
- chi -electron magnetic susceptibility tensor in cubic Angstroms
- ranges -Cartesian axis extents for the unpaired electron probability
- density cube as [xmin xmax ymin ymax zmin zmax] in Angstroms
- nxyz -nuclear coordinates as [x y z] with multiple rows) at which
- PCS is to be evaluated, in Angstroms
- Output:
- pcs_vals -pseudocontact shift in ppm at each nucleus
- pcs_cube -pseudocontact shift field on the same grid as the unpaired
