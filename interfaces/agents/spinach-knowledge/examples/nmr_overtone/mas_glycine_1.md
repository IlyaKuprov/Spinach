# examples/nmr_overtone/mas_glycine_1.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_overtone/mas_glycine_1.m`
- Signature: `mas_glycine_1()`
- Total lines: 76

## Purpose

Overtone detection 14N magic angle spinning NMR spectrum of glycine, computed using Fokker-Planck formalism. Glycine quadrupolar tensor data comes from the paper by O'Dell and Ratcliffe: Simulation parameters are set to reproduce Figure 3b from our paper on the subject: Calculation time: minutes

## Physical / mathematical content

- Overtone NMR examples. The important regime is excitation or detection of formally forbidden high-order transitions in quadrupolar nuclei, usually aided by MAS or Fokker-Planck treatments of periodic motion.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.
- Propagation is accelerated with a Krylov-subspace method, replacing direct matrix exponentiation by projection into a much smaller Arnoldi/Lanczos-type subspace.
- Quadrupolar physics is relevant: nuclei with spin > 1/2 interact with the electric field gradient tensor, introducing second-rank anisotropy, asymmetry, and overtone or MQ phenomena.

## Numerical / algorithmic content

- A Krylov-subspace or Arnoldi construction is used to avoid forming or exponentiating very large dense propagators directly.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 20-21: System specification; implemented by `sys.magnet=14.1; sys.isotopes={'14N'}`.
- Lines 25-26: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 29-30: Relaxation theory; implemented by `inter.relaxation={'damp'}`.
- Lines 35-36: Algorithmic options; implemented by `sys.disable={'krylov','trajlevel'}`.
- Lines 38-39: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 42-43: Magic angle; implemented by `theta=atan(sqrt(2))`.
- Lines 45-46: Spectrum setup; implemented by `parameters.max_rank=6`.
- Lines 66-67: Simulation; implemented by `spectrum=singlerot(spin_system,@overtone_pa,parameters,'qnmr')`.
- Lines 69-70: Phasing; implemented by `spectrum=exp(1i*1.35)*spectrum`.
- Lines 72-73: Plotting; implemented by `kfigure(); plot_1d(spin_system,real(spectrum),parameters)`.

### Key state/data transformations

- Lines 21: computes `sys.magnet` using `sys.magnet=14.1; sys.isotopes={'14N'}`.
- Lines 22: computes `inter.coupling.matrix{1,1}` using `inter.coupling.matrix{1,1}=eeqq2nqi(1.18e6,0.53,1,[0 0 0])`.
- Lines 23: computes `inter.zeeman.scalar{1}` using `inter.zeeman.scalar{1}=32.4`.
- Lines 26: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 27: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 30: computes `inter.relaxation` using `inter.relaxation={'damp'}`.
- Lines 31: computes `inter.rlx_keep` using `inter.rlx_keep='diagonal'`.
- Lines 32: computes `inter.equilibrium` using `inter.equilibrium='zero'`.
- Lines 33: computes `inter.damp_rate` using `inter.damp_rate=300`.
- Lines 36: computes `sys.disable` using `sys.disable={'krylov','trajlevel'}`.
- Lines 39: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 43: computes `theta` using `theta=atan(sqrt(2))`.
- Lines 46: computes `parameters.max_rank` using `parameters.max_rank=6`.
- Lines 47: computes `parameters.axis` using `parameters.axis=[sqrt(2/3) 0 sqrt(1/3)]`.
- Lines 48: computes `parameters.rate` using `parameters.rate=-19840`.
- Lines 49: computes `parameters.grid` using `parameters.grid='rep_2ang_6400pts_sph'`.
- Lines 50: computes `parameters.sweep` using `parameters.sweep=[44e3 52e3]`.
- Lines 51: computes `parameters.npoints` using `parameters.npoints=256`.

## Implementation structure

- Overtone detection 14N magic angle spinning NMR spectrum of glycine,
- computed using Fokker-Planck formalism. Glycine quadrupolar tensor
- data comes from the paper by O'Dell and Ratcliffe:
- Simulation parameters are set to reproduce Figure 3b from our paper
- on the subject:
- Calculation time: minutes
- System specification
- Basis set
- Relaxation theory
- Algorithmic options
- Spinach housekeeping
- Magic angle

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `eeqq2nqi()`, `create()`, `basis()`, `atan()`, `state()`, `operator()`, `singlerot()`, `kfigure()`, `plot_1d()`.
