# examples/nmr_overtone/dor_glycine.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_overtone/dor_glycine.m`
- Signature: `dor_glycine()`
- Total lines: 78

## Purpose

Panoramic double rotation overtone 14N spectrum of glycine, simulated as described in our paper (Figure 1B): A short pulse with instrumentally inaccessible power is gi- ven to make the excitation pattern uniform. Glycine quadru- polar tensor data comes from O'Dell and Ratcliffe: Calculation time: hours

## Physical / mathematical content

- Overtone NMR examples. The important regime is excitation or detection of formally forbidden high-order transitions in quadrupolar nuclei, usually aided by MAS or Fokker-Planck treatments of periodic motion.
- Propagation is accelerated with a Krylov-subspace method, replacing direct matrix exponentiation by projection into a much smaller Arnoldi/Lanczos-type subspace.

## Numerical / algorithmic content

- A Krylov-subspace or Arnoldi construction is used to avoid forming or exponentiating very large dense propagators directly.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 20-21: System specification; implemented by `sys.magnet=14.1; sys.isotopes={'14N'}`.
- Lines 24-25: Relaxation theory; implemented by `inter.relaxation={'damp'}`.
- Lines 30-31: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 34-35: Algorithmic options; implemented by `sys.disable={'krylov','trajlevel'}`.
- Lines 37-38: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 41-42: Magic angle; implemented by `theta=atan(sqrt(2))`.
- Lines 44-45: Experiment setup; implemented by `parameters.rate_outer=1425`.
- Lines 68-69: Run the simulation; implemented by `spectrum=doublerot(spin_system,@overtone_pa,parameters,'qnmr')`.
- Lines 71-72: Phasing; implemented by `spectrum=exp(1i*1.49)*spectrum`.
- Lines 74-75: Plotting; implemented by `kfigure(); plot_1d(spin_system,real(spectrum),parameters)`.

### Key state/data transformations

- Lines 21: computes `sys.magnet` using `sys.magnet=14.1; sys.isotopes={'14N'}`.
- Lines 22: computes `inter.coupling.matrix{1,1}` using `inter.coupling.matrix{1,1}=eeqq2nqi(1.18e6,0.53,1,[0 0 0])`.
- Lines 25: computes `inter.relaxation` using `inter.relaxation={'damp'}`.
- Lines 26: computes `inter.rlx_keep` using `inter.rlx_keep='diagonal'`.
- Lines 27: computes `inter.equilibrium` using `inter.equilibrium='zero'`.
- Lines 28: computes `inter.damp_rate` using `inter.damp_rate=100`.
- Lines 31: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 32: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 35: computes `sys.disable` using `sys.disable={'krylov','trajlevel'}`.
- Lines 38: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 42: computes `theta` using `theta=atan(sqrt(2))`.
- Lines 45: computes `parameters.rate_outer` using `parameters.rate_outer=1425`.
- Lines 46: computes `parameters.rate_inner` using `parameters.rate_inner=6950`.
- Lines 47: computes `parameters.rank_outer` using `parameters.rank_outer=5`.
- Lines 48: computes `parameters.rank_inner` using `parameters.rank_inner=5`.
- Lines 49: computes `parameters.axis_outer` using `parameters.axis_outer=[sin(theta) 0 cos(theta)]`.
- Lines 50: computes `parameters.axis_inner` using `parameters.axis_inner=[sqrt(20-2*sqrt(30)) 0 sqrt(15+2*sqrt(30))]`.
- Lines 51: computes `parameters.grid` using `parameters.grid='rep_2ang_200pts_sph'`.

## Implementation structure

- Panoramic double rotation overtone 14N spectrum of glycine,
- simulated as described in our paper (Figure 1B):
- A short pulse with instrumentally inaccessible power is gi-
- ven to make the excitation pattern uniform. Glycine quadru-
- polar tensor data comes from O'Dell and Ratcliffe:
- Calculation time: hours
- System specification
- Relaxation theory
- Basis set
- Algorithmic options
- Spinach housekeeping
- Magic angle

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `eeqq2nqi()`, `create()`, `basis()`, `atan()`, `state()`, `operator()`, `doublerot()`, `kfigure()`, `plot_1d()`.
