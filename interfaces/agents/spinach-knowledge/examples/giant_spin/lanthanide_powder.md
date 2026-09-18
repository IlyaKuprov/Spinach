# examples/giant_spin/lanthanide_powder.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/giant_spin/lanthanide_powder.m`
- Signature: `lanthanide_powder()`
- Total lines: 57

## Purpose

Powder spectrum of Gd(III) with ZFS up to 4th spherical rank using the giant spin Hamiltonian formalism in a sweepable 400 MHz NMR magnet and microwaves at 263.2 GHz. Odd ranks are zero because a zero-field Hamiltonian must be even under time reversal. The 4th rank terms are converted from the Stevens parameters b40=4e-4 cm^-1 and b44=-2e-4 cm^-1 reported for Gd(III) in tetragonal BaTiO3 by Rimai and deMars (https://doi.org/10.1103/PhysRev.127.702). Calculation time: seconds.

## Physical / mathematical content

- Giant-spin examples. The effective model treats lanthanides or high-spin centres using crystal-field / Stevens-operator Hamiltonians, Zeeman splitting, and magnetisation dynamics.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 17-18: Spin system properties; implemented by `sys.isotopes={'E8'}`.
- Lines 24-25: Field sweep; implemented by `sys.magnet=1`.
- Lines 27-28: Basis set; implemented by `bas.formalism='zeeman-hilb'`.
- Lines 31-32: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 35-36: Experiment parameters; implemented by `parameters.spins={'E8'}`.
- Lines 46-47: Run the simulation in the high-T approximation; implemented by `parameters.rho0=-state(spin_system,'Lz','E8')`.
- Lines 50-51: Plotting; implemented by `kfigure(); plot(parameters.b_axis,spec)`.

### Key state/data transformations

- Lines 18: computes `sys.isotopes` using `sys.isotopes={'E8'}`.
- Lines 19: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={1.9918}`.
- Lines 20-21: computes `inter.giant.coeff` using `inter.giant.coeff={{[0 0 0],[0 0 -4.65e8 0 0],[0 0 0 0 0 0 0],[-2.00e5 0 0 0 3.34e6 0 0 0 -2.00e5]}}`.
- Lines 22: computes `inter.giant.euler` using `inter.giant.euler={{[0 0 0],[0 0 0],[0 0 0],[0 0 0]}}`.
- Lines 25: computes `sys.magnet` using `sys.magnet=1`.
- Lines 28: computes `bas.formalism` using `bas.formalism='zeeman-hilb'`.
- Lines 29: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 32: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 36: computes `parameters.spins` using `parameters.spins={'E8'}`.
- Lines 37: computes `parameters.grid` using `parameters.grid='rep_2ang_100pts_sph'`.
- Lines 38: computes `parameters.mw_freq` using `parameters.mw_freq=263.2e9`.
- Lines 39: computes `parameters.fwhm` using `parameters.fwhm=2e-4`.
- Lines 40: computes `parameters.int_tol` using `parameters.int_tol=10.0`.
- Lines 41: computes `parameters.tm_tol` using `parameters.tm_tol=0.1`.
- Lines 42: computes `parameters.window` using `parameters.window=[9.32 9.56]`.
- Lines 43: computes `parameters.npoints` using `parameters.npoints=4096`.
- Lines 44: computes `parameters.rspt_order` using `parameters.rspt_order=Inf`.
- Lines 47: computes `parameters.rho0` using `parameters.rho0=-state(spin_system,'Lz','E8')`.

## Implementation structure

- Powder spectrum of Gd(III) with ZFS up to 4th spherical rank
- using the giant spin Hamiltonian formalism in a sweepable 400
- MHz NMR magnet and microwaves at 263.2 GHz. Odd ranks are zero
- because a zero-field Hamiltonian must be even under time rever-
- sal. The 4th rank terms are converted from the Stevens parame-
- ters b40=4e-4 cm^-1 and b44=-2e-4 cm^-1 reported for Gd(III)
- in tetragonal BaTiO3 by Rimai and deMars:
- Calculation time: seconds.
- Spin system properties
- Field sweep
- Basis set
- Spinach housekeeping
- Experiment parameters
- Run the simulation in the high-T approximation
- Plotting

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `fieldsweep()`, `kfigure()`, `kxlabel()`, `kylabel()`.
