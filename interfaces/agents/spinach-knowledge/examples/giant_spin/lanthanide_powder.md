# examples/giant_spin/lanthanide_powder.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/giant_spin/lanthanide_powder.m`
- Signature: `lanthanide_powder()`
- Total lines: 50

## Purpose

Powder spectrum of Gd(III) with ZFS up to 3rd spherical rank using the giant spin Hamiltonian formalism in a sweepable 400 MHz NMR magnet and microwaves at 263.2 GHz. Calculation time: seconds.

## Physical / mathematical content

- Giant-spin examples. The effective model treats lanthanides or high-spin centres using crystal-field / Stevens-operator Hamiltonians, Zeeman splitting, and magnetisation dynamics.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 11-12: Spin system properties; implemented by `sys.isotopes={'E8'}`.
- Lines 17-18: Field sweep; implemented by `sys.magnet=1`.
- Lines 20-21: Basis set; implemented by `bas.formalism='zeeman-hilb'`.
- Lines 24-25: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 28-29: Experiment parameters; implemented by `parameters.spins={'E8'}`.
- Lines 39-40: Run the simulation in the high-T approximation; implemented by `parameters.rho0=-state(spin_system,'Lz','E8')`.
- Lines 43-44: Plotting; implemented by `kfigure(); plot(parameters.b_axis,spec)`.

### Key state/data transformations

- Lines 12: computes `sys.isotopes` using `sys.isotopes={'E8'}`.
- Lines 13: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={1.9918}`.
- Lines 14: computes `inter.giant.coeff` using `inter.giant.coeff={{[0 0 0],[0 0 -4.65e8 0 0],[1e7 0 0 2e7 0 0 1e7]}}`.
- Lines 15: computes `inter.giant.euler` using `inter.giant.euler={{[0 0 0],[0 0 0],[0 0 0]}}`.
- Lines 18: computes `sys.magnet` using `sys.magnet=1`.
- Lines 21: computes `bas.formalism` using `bas.formalism='zeeman-hilb'`.
- Lines 22: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 25: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 29: computes `parameters.spins` using `parameters.spins={'E8'}`.
- Lines 30: computes `parameters.grid` using `parameters.grid='rep_2ang_100pts_sph'`.
- Lines 31: computes `parameters.mw_freq` using `parameters.mw_freq=263.2e9`.
- Lines 32: computes `parameters.fwhm` using `parameters.fwhm=2e-4`.
- Lines 33: computes `parameters.int_tol` using `parameters.int_tol=10.0`.
- Lines 34: computes `parameters.tm_tol` using `parameters.tm_tol=0.1`.
- Lines 35: computes `parameters.window` using `parameters.window=[9.32 9.56]`.
- Lines 36: computes `parameters.npoints` using `parameters.npoints=4096`.
- Lines 37: computes `parameters.rspt_order` using `parameters.rspt_order=Inf`.
- Lines 40: computes `parameters.rho0` using `parameters.rho0=-state(spin_system,'Lz','E8')`.

## Implementation structure

- Powder spectrum of Gd(III) with ZFS up to 3rd spherical rank
- using the giant spin Hamiltonian formalism in a sweepable 400
- MHz NMR magnet and microwaves at 263.2 GHz.
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
