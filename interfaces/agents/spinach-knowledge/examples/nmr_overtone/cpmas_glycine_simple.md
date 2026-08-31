# examples/nmr_overtone/cpmas_glycine_simple.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_overtone/cpmas_glycine_simple.m`
- Signature: `cpmas_glycine_simple()`
- Total lines: 74

## Purpose

Cross-polarization experiment between protons and 14N overtone transition in glycine under MAS. Glycine quadrupolar tensor da- ta comes from the paper by O'Dell and Ratcliffe: Calculation time: hours

## Physical / mathematical content

- Overtone NMR examples. The important regime is excitation or detection of formally forbidden high-order transitions in quadrupolar nuclei, usually aided by MAS or Fokker-Planck treatments of periodic motion.
- Propagation is accelerated with a Krylov-subspace method, replacing direct matrix exponentiation by projection into a much smaller Arnoldi/Lanczos-type subspace.
- Quadrupolar physics is relevant: nuclei with spin > 1/2 interact with the electric field gradient tensor, introducing second-rank anisotropy, asymmetry, and overtone or MQ phenomena.

## Numerical / algorithmic content

- A Krylov-subspace or Arnoldi construction is used to avoid forming or exponentiating very large dense propagators directly.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 15-16: System specification; implemented by `sys.magnet=14.1; sys.isotopes={'14N','1H'}`.
- Lines 23-24: Relaxation theory; implemented by `inter.relaxation={'damp'}`.
- Lines 29-30: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 33-34: Algorithmic options; implemented by `sys.disable={'krylov','trajlevel'}`.
- Lines 36-37: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 40-41: Magic angle; implemented by `theta=atan(sqrt(2))`.
- Lines 43-44: Spectrum setup; implemented by `parameters.max_rank=7`.
- Lines 67-68: Simulation; implemented by `spectrum=singlerot(spin_system,@overtone_cp,parameters,'qnmr')`.
- Lines 70-71: Plotting; implemented by `kfigure(); plot_1d(spin_system,real(spectrum),parameters)`.

### Key state/data transformations

- Lines 16: computes `sys.magnet` using `sys.magnet=14.1; sys.isotopes={'14N','1H'}`.
- Lines 17: computes `inter.coupling.matrix{1,1}` using `inter.coupling.matrix{1,1}=eeqq2nqi(1.18e6,0.53,1,[0 0 0])`.
- Lines 18: computes `inter.coupling.matrix{2,2}` using `inter.coupling.matrix{2,2}=[]`.
- Lines 19: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={32.4 0}`.
- Lines 20: computes `inter.coordinates` using `inter.coordinates={[0 0 0]`.
- Lines 24: computes `inter.relaxation` using `inter.relaxation={'damp'}`.
- Lines 25: computes `inter.rlx_keep` using `inter.rlx_keep='diagonal'`.
- Lines 26: computes `inter.equilibrium` using `inter.equilibrium='zero'`.
- Lines 27: computes `inter.damp_rate` using `inter.damp_rate=300`.
- Lines 30: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 31: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 34: computes `sys.disable` using `sys.disable={'krylov','trajlevel'}`.
- Lines 37: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 41: computes `theta` using `theta=atan(sqrt(2))`.
- Lines 44: computes `parameters.max_rank` using `parameters.max_rank=7`.
- Lines 45: computes `parameters.axis` using `parameters.axis=[sqrt(2/3) 0 sqrt(1/3)]`.
- Lines 46: computes `parameters.rate` using `parameters.rate=-19840`.
- Lines 47: computes `parameters.grid` using `parameters.grid='rep_2ang_6400pts_sph'`.

## Implementation structure

- Cross-polarization experiment between protons and 14N overtone
- transition in glycine under MAS. Glycine quadrupolar tensor da-
- ta comes from the paper by O'Dell and Ratcliffe:
- Calculation time: hours
- System specification
- Relaxation theory
- Basis set
- Algorithmic options
- Spinach housekeeping
- Magic angle
- Spectrum setup
- Simulation

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `eeqq2nqi()`, `create()`, `basis()`, `atan()`, `state()`, `operator()`, `singlerot()`, `kfigure()`, `plot_1d()`.
