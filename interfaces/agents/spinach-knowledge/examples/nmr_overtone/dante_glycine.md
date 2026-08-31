# examples/nmr_overtone/dante_glycine.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_overtone/dante_glycine.m`
- Signature: `dante_glycine()`
- Total lines: 73

## Purpose

14N overtone DANTE spectrum of glycine, computed using Fokker-Planck formalism. Glycine quadrupolar tensor data comes from the paper by O'Dell and Ratcliffe: Calculation time: minutes

## Physical / mathematical content

- Overtone NMR examples. The important regime is excitation or detection of formally forbidden high-order transitions in quadrupolar nuclei, usually aided by MAS or Fokker-Planck treatments of periodic motion.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.
- Propagation is accelerated with a Krylov-subspace method, replacing direct matrix exponentiation by projection into a much smaller Arnoldi/Lanczos-type subspace.
- Quadrupolar physics is relevant: nuclei with spin > 1/2 interact with the electric field gradient tensor, introducing second-rank anisotropy, asymmetry, and overtone or MQ phenomena.

## Numerical / algorithmic content

- A Krylov-subspace or Arnoldi construction is used to avoid forming or exponentiating very large dense propagators directly.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 15-16: System specification; implemented by `sys.magnet=14.1; sys.isotopes={'14N'}`.
- Lines 20-21: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 24-25: Relaxation theory; implemented by `inter.relaxation={'damp'}`.
- Lines 30-31: Algorithmic options; implemented by `sys.disable={'krylov','trajlevel'}`.
- Lines 33-34: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 37-38: Magic angle; implemented by `theta=atan(sqrt(2))`.
- Lines 40-41: Spectrum setup; implemented by `parameters.max_rank=7`.
- Lines 63-64: Simulation; implemented by `spectrum=singlerot(spin_system,@overtone_dante,parameters,'qnmr')`.
- Lines 66-67: Phasing; implemented by `spectrum=exp(-1i*2.12)*spectrum`.
- Lines 69-70: Plotting; implemented by `kfigure(); plot_1d(spin_system,real(spectrum),parameters)`.

### Key state/data transformations

- Lines 16: computes `sys.magnet` using `sys.magnet=14.1; sys.isotopes={'14N'}`.
- Lines 17: computes `inter.coupling.matrix{1,1}` using `inter.coupling.matrix{1,1}=eeqq2nqi(1.18e6,0.53,1,[0 0 0])`.
- Lines 18: computes `inter.zeeman.scalar{1}` using `inter.zeeman.scalar{1}=32.4`.
- Lines 21: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 22: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 25: computes `inter.relaxation` using `inter.relaxation={'damp'}`.
- Lines 26: computes `inter.rlx_keep` using `inter.rlx_keep='diagonal'`.
- Lines 27: computes `inter.equilibrium` using `inter.equilibrium='zero'`.
- Lines 28: computes `inter.damp_rate` using `inter.damp_rate=300`.
- Lines 31: computes `sys.disable` using `sys.disable={'krylov','trajlevel'}`.
- Lines 34: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 38: computes `theta` using `theta=atan(sqrt(2))`.
- Lines 41: computes `parameters.max_rank` using `parameters.max_rank=7`.
- Lines 42: computes `parameters.axis` using `parameters.axis=[sqrt(2/3) 0 sqrt(1/3)]`.
- Lines 43: computes `parameters.rate` using `parameters.rate=-19840`.
- Lines 44: computes `parameters.grid` using `parameters.grid='rep_2ang_1600pts_sph'`.
- Lines 45: computes `parameters.sweep` using `parameters.sweep=[-60e3 80e3]`.
- Lines 46: computes `parameters.npoints` using `parameters.npoints=2048`.

## Implementation structure

- 14N overtone DANTE spectrum of glycine, computed using Fokker-Planck
- formalism. Glycine quadrupolar tensor data comes from the paper by
- O'Dell and Ratcliffe:
- Calculation time: minutes
- System specification
- Basis set
- Relaxation theory
- Algorithmic options
- Spinach housekeeping
- Magic angle
- Spectrum setup
- Simulation

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `eeqq2nqi()`, `create()`, `basis()`, `atan()`, `state()`, `operator()`, `singlerot()`, `kfigure()`, `plot_1d()`.
