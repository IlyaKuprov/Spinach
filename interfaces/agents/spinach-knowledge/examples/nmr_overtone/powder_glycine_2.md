# examples/nmr_overtone/powder_glycine_2.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_overtone/powder_glycine_2.m`
- Signature: `powder_glycine_2()`
- Total lines: 60

## Purpose

Overtone detection 14N powder NMR spectrum of glycine, computed using Fokker-Planck formalism. Glycine quadrupolar tensor data comes from the paper by O'Dell and Ratcliffe: This simulation demonstrates that the spin state that gives rise to the overtone signal in a static sample is the T2,-2 coherence. Calculation time: seconds

## Physical / mathematical content

- Overtone NMR examples. The important regime is excitation or detection of formally forbidden high-order transitions in quadrupolar nuclei, usually aided by MAS or Fokker-Planck treatments of periodic motion.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.
- Propagation is accelerated with a Krylov-subspace method, replacing direct matrix exponentiation by projection into a much smaller Arnoldi/Lanczos-type subspace.
- Quadrupolar physics is relevant: nuclei with spin > 1/2 interact with the electric field gradient tensor, introducing second-rank anisotropy, asymmetry, and overtone or MQ phenomena.

## Numerical / algorithmic content

- A Krylov-subspace or Arnoldi construction is used to avoid forming or exponentiating very large dense propagators directly.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 16-17: System specification; implemented by `sys.magnet=14.1; sys.isotopes={'14N'}`.
- Lines 21-22: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 25-26: Relaxation theory; implemented by `inter.relaxation={'damp'}`.
- Lines 31-32: Algorithmic options; implemented by `sys.disable={'krylov','trajlevel'}`.
- Lines 34-35: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 38-39: Magic angle; implemented by `theta=atan(sqrt(2))`.
- Lines 41-42: Spectrum setup; implemented by `parameters.grid='rep_2ang_6400pts_sph'`.
- Lines 53-54: Simulation; implemented by `spectrum=powder(spin_system,@overtone_a,parameters,'qnmr')`.
- Lines 56-57: Plotting; implemented by `kfigure(); plot_1d(spin_system,real(spectrum),parameters)`.

### Key state/data transformations

- Lines 17: computes `sys.magnet` using `sys.magnet=14.1; sys.isotopes={'14N'}`.
- Lines 18: computes `inter.coupling.matrix{1,1}` using `inter.coupling.matrix{1,1}=eeqq2nqi(1.18e6,0.53,1,[0 0 0])`.
- Lines 19: computes `inter.zeeman.scalar{1}` using `inter.zeeman.scalar{1}=32.4`.
- Lines 22: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 23: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 26: computes `inter.relaxation` using `inter.relaxation={'damp'}`.
- Lines 27: computes `inter.rlx_keep` using `inter.rlx_keep='diagonal'`.
- Lines 28: computes `inter.equilibrium` using `inter.equilibrium='zero'`.
- Lines 29: computes `inter.damp_rate` using `inter.damp_rate=500`.
- Lines 32: computes `sys.disable` using `sys.disable={'krylov','trajlevel'}`.
- Lines 35: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 39: computes `theta` using `theta=atan(sqrt(2))`.
- Lines 42: computes `parameters.grid` using `parameters.grid='rep_2ang_6400pts_sph'`.
- Lines 43: computes `parameters.sweep` using `parameters.sweep=[0e3 15e3]`.
- Lines 44: computes `parameters.npoints` using `parameters.npoints=256`.
- Lines 45: computes `parameters.zerofill` using `parameters.zerofill=256`.
- Lines 46: computes `parameters.rho0` using `parameters.rho0=state(spin_system,'T2,-2','14N')`.
- Lines 47-48: computes `parameters.coil` using `parameters.coil=cos(theta)*state(spin_system,'Lz','14N')+ sin(theta)*state(spin_system,'Lx','14N')`.

## Implementation structure

- Overtone detection 14N powder NMR spectrum of glycine, computed using
- Fokker-Planck formalism. Glycine quadrupolar tensor data comes from
- the paper by O'Dell and Ratcliffe:
- This simulation demonstrates that the spin state that gives rise to
- the overtone signal in a static sample is the T2,-2 coherence.
- Calculation time: seconds
- System specification
- Basis set
- Relaxation theory
- Algorithmic options
- Spinach housekeeping
- Magic angle

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `eeqq2nqi()`, `create()`, `basis()`, `atan()`, `state()`, `powder()`, `kfigure()`, `plot_1d()`.
