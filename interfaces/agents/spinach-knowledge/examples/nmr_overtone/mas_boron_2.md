# examples/nmr_overtone/mas_boron_2.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_overtone/mas_boron_2.m`
- Signature: `mas_boron_2()`
- Total lines: 66

## Purpose

Overtone 10B magic angle spinning NMR spectrum. The sample is spinning in the JEOL direction. Parameters from Nghia Duong and Yusuke Nishiyama, realistic RF power and pulse width. Calculation time: hours

## Physical / mathematical content

- Overtone NMR examples. The important regime is excitation or detection of formally forbidden high-order transitions in quadrupolar nuclei, usually aided by MAS or Fokker-Planck treatments of periodic motion.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 11-12: System specification; implemented by `sys.magnet=16.4; sys.isotopes={'10B'}`.
- Lines 15-16: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 19-20: Relaxation theory; implemented by `inter.relaxation={'damp'}`.
- Lines 25-26: Algorithmic options; implemented by `sys.disable={'trajlevel'}`.
- Lines 28-29: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 32-33: Magic angle; implemented by `theta=atan(sqrt(2))`.
- Lines 35-36: Sequence parameters; implemented by `parameters.max_rank=12`.
- Lines 56-57: Simulation; implemented by `spectrum=singlerot(spin_system,@overtone_pa,parameters,'qnmr')`.
- Lines 59-60: Phasing; implemented by `spectrum=exp(1i*1.45)*spectrum`.
- Lines 62-63: Plotting; implemented by `kfigure(); plot_1d(spin_system,real(spectrum),parameters)`.

### Key state/data transformations

- Lines 12: computes `sys.magnet` using `sys.magnet=16.4; sys.isotopes={'10B'}`.
- Lines 13: computes `inter.coupling.matrix{1,1}` using `inter.coupling.matrix{1,1}=eeqq2nqi(0.7e6,0.0,3,[0 0 0])`.
- Lines 16: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 17: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 20: computes `inter.relaxation` using `inter.relaxation={'damp'}`.
- Lines 21: computes `inter.rlx_keep` using `inter.rlx_keep='diagonal'`.
- Lines 22: computes `inter.equilibrium` using `inter.equilibrium='zero'`.
- Lines 23: computes `inter.damp_rate` using `inter.damp_rate=100`.
- Lines 26: computes `sys.disable` using `sys.disable={'trajlevel'}`.
- Lines 29: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 33: computes `theta` using `theta=atan(sqrt(2))`.
- Lines 36: computes `parameters.max_rank` using `parameters.max_rank=12`.
- Lines 37: computes `parameters.axis` using `parameters.axis=[sqrt(2/3) 0 sqrt(1/3)]`.
- Lines 38: computes `parameters.rate` using `parameters.rate=70000`.
- Lines 39: computes `parameters.grid` using `parameters.grid='rep_2ang_200pts_sph'`.
- Lines 40: computes `parameters.sweep` using `parameters.sweep=[-141e3 -139e3]`.
- Lines 41: computes `parameters.npoints` using `parameters.npoints=256`.
- Lines 42: computes `parameters.zerofill` using `parameters.zerofill=256`.

## Implementation structure

- Overtone 10B magic angle spinning NMR spectrum. The sample is
- spinning in the JEOL direction. Parameters from Nghia Duong
- and Yusuke Nishiyama, realistic RF power and pulse width.
- Calculation time: hours
- System specification
- Basis set
- Relaxation theory
- Algorithmic options
- Spinach housekeeping
- Magic angle
- Sequence parameters
- Simulation

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `eeqq2nqi()`, `create()`, `basis()`, `atan()`, `state()`, `operator()`, `singlerot()`, `kfigure()`, `plot_1d()`.
