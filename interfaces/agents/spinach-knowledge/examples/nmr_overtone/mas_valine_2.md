# examples/nmr_overtone/mas_valine_2.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_overtone/mas_valine_2.m`
- Signature: `mas_valine_2()`
- Total lines: 59

## Purpose

Overtone Z-detection 14N magic angle spinning NMR spectrum of N-acetylvaline, computed using Fokker-Planck formalism. Valine quadrupolar tensor data comes from our paper: Calculation time: hours

## Physical / mathematical content

- Overtone NMR examples. The important regime is excitation or detection of formally forbidden high-order transitions in quadrupolar nuclei, usually aided by MAS or Fokker-Planck treatments of periodic motion.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.
- Propagation is accelerated with a Krylov-subspace method, replacing direct matrix exponentiation by projection into a much smaller Arnoldi/Lanczos-type subspace.
- Quadrupolar physics is relevant: nuclei with spin > 1/2 interact with the electric field gradient tensor, introducing second-rank anisotropy, asymmetry, and overtone or MQ phenomena.

## Numerical / algorithmic content

- A Krylov-subspace or Arnoldi construction is used to avoid forming or exponentiating very large dense propagators directly.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 15-16: System specification; implemented by `sys.magnet=14.102; sys.isotopes={'14N'}`.
- Lines 21-22: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 25-26: Relaxation theory; implemented by `inter.relaxation={'damp'}`.
- Lines 31-32: Algorithmic options; implemented by `sys.disable={'krylov','trajlevel'}`.
- Lines 34-35: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 38-39: Spectrum setup; implemented by `parameters.max_rank=8`.
- Lines 52-53: Simulation; implemented by `spectrum=singlerot(spin_system,@overtone_a,parameters,'qnmr')`.
- Lines 55-56: Plotting; implemented by `kfigure(); plot_1d(spin_system,real(spectrum),parameters)`.

### Key state/data transformations

- Lines 16: computes `sys.magnet` using `sys.magnet=14.102; sys.isotopes={'14N'}`.
- Lines 17: computes `inter.coupling.matrix{1,1}` using `inter.coupling.matrix{1,1}=eeqq2nqi(3.21e6,0.27,1,[0 0 0])`.
- Lines 18: computes `inter.zeeman.eigs` using `inter.zeeman.eigs={[57.5 81.0 227.0]}`.
- Lines 19: computes `inter.zeeman.euler` using `inter.zeeman.euler={[-90 -90 -17]*(pi/180)}`.
- Lines 22: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 23: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 26: computes `inter.relaxation` using `inter.relaxation={'damp'}`.
- Lines 27: computes `inter.rlx_keep` using `inter.rlx_keep='diagonal'`.
- Lines 28: computes `inter.equilibrium` using `inter.equilibrium='zero'`.
- Lines 29: computes `inter.damp_rate` using `inter.damp_rate=2000`.
- Lines 32: computes `sys.disable` using `sys.disable={'krylov','trajlevel'}`.
- Lines 35: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 39: computes `parameters.max_rank` using `parameters.max_rank=8`.
- Lines 40: computes `parameters.axis` using `parameters.axis=[sqrt(2/3) 0 sqrt(1/3)]`.
- Lines 41: computes `parameters.rate` using `parameters.rate=-19840`.
- Lines 42: computes `parameters.grid` using `parameters.grid='rep_2ang_6400pts_sph'`.
- Lines 43: computes `parameters.sweep` using `parameters.sweep=[75e3 100e3]`.
- Lines 44: computes `parameters.npoints` using `parameters.npoints=256`.

## Implementation structure

- Overtone Z-detection 14N magic angle spinning NMR spectrum
- of N-acetylvaline, computed using Fokker-Planck formalism.
- Valine quadrupolar tensor data comes from our paper:
- Calculation time: hours
- System specification
- Basis set
- Relaxation theory
- Algorithmic options
- Spinach housekeeping
- Spectrum setup
- Simulation
- Plotting

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `eeqq2nqi()`, `create()`, `basis()`, `state()`, `singlerot()`, `kfigure()`, `plot_1d()`.
