# examples/nmr_solids/redor_curve.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_solids/redor_curve.m`
- Signature: `redor_curve()`
- Total lines: 54

## Purpose

REDOR dephasing curve for a simple 13C-15N spin pair using Fokker-Planck MAS formalism. Calculation time: seconds

## Physical / mathematical content

- Solid-state NMR examples. The key physics is anisotropic spin interactions under static or magic-angle-spinning conditions: chemical-shift anisotropy, dipolar coupling, quadrupolar coupling, cross-polarisation, and orientation averaging using Floquet, Fokker-Planck, or direct powder quadrature formalisms.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 10-11: System specification; implemented by `sys.magnet=14.1`.
- Lines 14-15: Interactions; implemented by `inter.zeeman.scalar={0.0 0.0}`.
- Lines 19-20: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 23-24: Algorithmic options; implemented by `sys.disable={'trajlevel'}`.
- Lines 27-28: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 31-32: REDOR setup; implemented by `parameters.spins={'13C','15N'}`.
- Lines 42-43: Simulation; implemented by `curve=singlerot(spin_system,@redor,parameters,'nmr')`.
- Lines 45-46: Normalised REDOR difference; implemented by `redor_diff=real(curve(3,:)./curve(1,:))`.
- Lines 48-49: Plotting; implemented by `kfigure(); plot(parameters.ncycles,redor_diff)`.

### Key state/data transformations

- Lines 11: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 12: computes `sys.isotopes` using `sys.isotopes={'13C','15N'}`.
- Lines 15: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={0.0 0.0}`.
- Lines 16: computes `inter.coordinates` using `inter.coordinates={[0.00 0.00 0.00]`.
- Lines 20: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 21: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 24: computes `sys.disable` using `sys.disable={'trajlevel'}`.
- Lines 25: computes `sys.enable` using `sys.enable={'prop_cache'}`.
- Lines 28: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 32: computes `parameters.spins` using `parameters.spins={'13C','15N'}`.
- Lines 33: computes `parameters.rate` using `parameters.rate=10000`.
- Lines 34: computes `parameters.axis` using `parameters.axis=[sqrt(2/3) 0 sqrt(1/3)]`.
- Lines 35: computes `parameters.max_rank` using `parameters.max_rank=9`.
- Lines 36: computes `parameters.grid` using `parameters.grid='leb_2ang_rank_23'`.
- Lines 37: computes `parameters.ncycles` using `parameters.ncycles=0:48`.
- Lines 38: computes `parameters.rho0` using `parameters.rho0=state(spin_system,'Lx','13C')`.
- Lines 39: computes `parameters.coil` using `parameters.coil=state(spin_system,'Lx','13C')`.
- Lines 40: computes `parameters.verbose` using `parameters.verbose=0`.

## Implementation structure

- REDOR dephasing curve for a simple 13C-15N spin pair using
- Fokker-Planck MAS formalism.
- Calculation time: seconds
- System specification
- Interactions
- Basis set
- Algorithmic options
- Spinach housekeeping
- REDOR setup
- Simulation
- Normalised REDOR difference
- Plotting

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `singlerot()`, `curve()`, `kfigure()`, `kxlabel()`, `kylabel()`.
