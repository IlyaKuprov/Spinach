# examples/nmr_overtone/mas_boron_1.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_overtone/mas_boron_1.m`
- Signature: `mas_boron_1()`
- Total lines: 54

## Purpose

Overtone Z-detection 10B magic angle spinning NMR spectrum. The sample is spinning in the JEOL direction. Parameters from Nghia Duong and Yusuke Nishiyama. The simulation focuses on the most intense of the five overtone spinning sidebands. Calculation time: hours

## Physical / mathematical content

- Overtone NMR examples. The important regime is excitation or detection of formally forbidden high-order transitions in quadrupolar nuclei, usually aided by MAS or Fokker-Planck treatments of periodic motion.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 12-13: System specification; implemented by `sys.magnet=16.4; sys.isotopes={'10B'}`.
- Lines 16-17: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 20-21: Relaxation theory; implemented by `inter.relaxation={'damp'}`.
- Lines 26-27: Algorithmic options; implemented by `sys.disable={'trajlevel'}`.
- Lines 29-30: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 33-34: Sequence parameters; implemented by `parameters.max_rank=12`.
- Lines 47-48: Simulation; implemented by `spectrum=singlerot(spin_system,@overtone_a,parameters,'qnmr')`.
- Lines 50-51: Plotting; implemented by `kfigure(); plot_1d(spin_system,real(spectrum),parameters)`.

### Key state/data transformations

- Lines 13: computes `sys.magnet` using `sys.magnet=16.4; sys.isotopes={'10B'}`.
- Lines 14: computes `inter.coupling.matrix{1,1}` using `inter.coupling.matrix{1,1}=eeqq2nqi(0.7e6,0.0,3,[0 0 0])`.
- Lines 17: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 18: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 21: computes `inter.relaxation` using `inter.relaxation={'damp'}`.
- Lines 22: computes `inter.rlx_keep` using `inter.rlx_keep='diagonal'`.
- Lines 23: computes `inter.equilibrium` using `inter.equilibrium='zero'`.
- Lines 24: computes `inter.damp_rate` using `inter.damp_rate=50`.
- Lines 27: computes `sys.disable` using `sys.disable={'trajlevel'}`.
- Lines 30: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 34: computes `parameters.max_rank` using `parameters.max_rank=12`.
- Lines 35: computes `parameters.axis` using `parameters.axis=[sqrt(2/3) 0 sqrt(1/3)]`.
- Lines 36: computes `parameters.rate` using `parameters.rate=70000`.
- Lines 37: computes `parameters.grid` using `parameters.grid='rep_2ang_800pts_sph'`.
- Lines 38: computes `parameters.sweep` using `parameters.sweep=[-141e3 -139e3]`.
- Lines 39: computes `parameters.npoints` using `parameters.npoints=256`.
- Lines 40: computes `parameters.zerofill` using `parameters.zerofill=256`.
- Lines 41: computes `parameters.rho0` using `parameters.rho0=state(spin_system,'Lz','10B')`.

## Implementation structure

- Overtone Z-detection 10B magic angle spinning NMR spectrum.
- The sample is spinning in the JEOL direction. Parameters from
- Nghia Duong and Yusuke Nishiyama. The simulation focuses on
- the most intense of the five overtone spinning sidebands.
- Calculation time: hours
- System specification
- Basis set
- Relaxation theory
- Algorithmic options
- Spinach housekeeping
- Sequence parameters
- Simulation

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `eeqq2nqi()`, `create()`, `basis()`, `state()`, `singlerot()`, `kfigure()`, `plot_1d()`.
