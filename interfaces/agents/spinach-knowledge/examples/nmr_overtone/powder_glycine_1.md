# examples/nmr_overtone/powder_glycine_1.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_overtone/powder_glycine_1.m`
- Signature: `powder_glycine_1()`
- Total lines: 65

## Purpose

Overtone detection 14N powder NMR spectrum of glycine, computed using Fokker-Planck formalism. Glycine quadrupolar tensor data comes from the paper by O'Dell and Ratcliffe: A very short pulse with an unphysically large power is used. Calculation time: seconds

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
- Lines 40-41: Spectrum setup; implemented by `parameters.grid='rep_2ang_6400pts_sph'`.
- Lines 58-59: Simulation; implemented by `spectrum=powder(spin_system,@overtone_pa,parameters,'qnmr')`.
- Lines 61-62: Plotting; implemented by `kfigure(); plot_1d(spin_system,real(spectrum),parameters)`.

### Key state/data transformations

- Lines 16: computes `sys.magnet` using `sys.magnet=14.1; sys.isotopes={'14N'}`.
- Lines 17: computes `inter.coupling.matrix{1,1}` using `inter.coupling.matrix{1,1}=eeqq2nqi(1.18e6,0.53,1,[0 0 0])`.
- Lines 18: computes `inter.zeeman.scalar{1}` using `inter.zeeman.scalar{1}=32.4`.
- Lines 21: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 22: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 25: computes `inter.relaxation` using `inter.relaxation={'damp'}`.
- Lines 26: computes `inter.rlx_keep` using `inter.rlx_keep='diagonal'`.
- Lines 27: computes `inter.equilibrium` using `inter.equilibrium='zero'`.
- Lines 28: computes `inter.damp_rate` using `inter.damp_rate=500`.
- Lines 31: computes `sys.disable` using `sys.disable={'krylov','trajlevel'}`.
- Lines 34: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 38: computes `theta` using `theta=atan(sqrt(2))`.
- Lines 41: computes `parameters.grid` using `parameters.grid='rep_2ang_6400pts_sph'`.
- Lines 42: computes `parameters.sweep` using `parameters.sweep=[0e3 15e3]`.
- Lines 43: computes `parameters.npoints` using `parameters.npoints=256`.
- Lines 44: computes `parameters.zerofill` using `parameters.zerofill=256`.
- Lines 45: computes `parameters.rho0` using `parameters.rho0=state(spin_system,'Lz','14N')`.
- Lines 46-47: computes `parameters.coil` using `parameters.coil=cos(theta)*state(spin_system,'Lz','14N')+ sin(theta)*state(spin_system,'Lx','14N')`.

## Implementation structure

- Overtone detection 14N powder NMR spectrum of glycine, computed using
- Fokker-Planck formalism. Glycine quadrupolar tensor data comes from
- the paper by O'Dell and Ratcliffe:
- A very short pulse with an unphysically large power is used.
- Calculation time: seconds
- System specification
- Basis set
- Relaxation theory
- Algorithmic options
- Spinach housekeeping
- Magic angle
- Spectrum setup

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `eeqq2nqi()`, `create()`, `basis()`, `atan()`, `state()`, `operator()`, `powder()`, `kfigure()`, `plot_1d()`.
