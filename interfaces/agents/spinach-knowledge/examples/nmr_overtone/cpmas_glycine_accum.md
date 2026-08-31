# examples/nmr_overtone/cpmas_glycine_accum.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_overtone/cpmas_glycine_accum.m`
- Signature: `cpmas_glycine_accum()`
- Total lines: 93

## Purpose

Cross-polarization experiment between protons and 14N overtone transition in glycine under MAS. Glycine quadrupolar tensor da- ta comes from the paper by O'Dell and Ratcliffe: Magnetisation accumulation profile as a function of contact ti- me with a rough powder grid. Calculation time: minutes

## Physical / mathematical content

- Overtone NMR examples. The important regime is excitation or detection of formally forbidden high-order transitions in quadrupolar nuclei, usually aided by MAS or Fokker-Planck treatments of periodic motion.
- Propagation is accelerated with a Krylov-subspace method, replacing direct matrix exponentiation by projection into a much smaller Arnoldi/Lanczos-type subspace.
- Quadrupolar physics is relevant: nuclei with spin > 1/2 interact with the electric field gradient tensor, introducing second-rank anisotropy, asymmetry, and overtone or MQ phenomena.

## Numerical / algorithmic content

- A Krylov-subspace or Arnoldi construction is used to avoid forming or exponentiating very large dense propagators directly.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 18-19: System specification; implemented by `sys.magnet=14.10220742; sys.isotopes={'14N','1H'}`.
- Lines 26-27: Relaxation theory; implemented by `inter.relaxation={'damp'}`.
- Lines 32-33: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 36-37: Algorithmic options; implemented by `sys.disable={'krylov','trajlevel'}`.
- Lines 39-40: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 43-44: Magic angle; implemented by `theta=atan(sqrt(2))`.
- Lines 46-47: Spectrum setup; implemented by `parameters.max_rank=7`.
- Lines 69-70: Start a new figure; implemented by `kfigure(); scale_figure([2.5 1.0])`.
- Lines 72-73: Loop over contact times; implemented by `for n=1:10`.
- Lines 75-76: Subplot selection; implemented by `subplot(1,10,n)`.
- Lines 78-79: Contact time; implemented by `parameters.rf_dur=1e-5*n`.
- Lines 81-82: Simulation; implemented by `spectrum=singlerot(spin_system,@overtone_cp,parameters,'qnmr')`.
- Lines 84-85: Plotting; implemented by `plot_1d(spin_system,real(spectrum),parameters)`.

### Control flow inferred from the code

- Line 73: `for` loop over `n=1:10`.

### Key state/data transformations

- Lines 19: computes `sys.magnet` using `sys.magnet=14.10220742; sys.isotopes={'14N','1H'}`.
- Lines 20: computes `inter.coupling.matrix{1,1}` using `inter.coupling.matrix{1,1}=eeqq2nqi(1.18e6,0.53,1,[0 0 0])`.
- Lines 21: computes `inter.coupling.matrix{2,2}` using `inter.coupling.matrix{2,2}=[]`.
- Lines 22: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={32.4 0}`.
- Lines 23: computes `inter.coordinates` using `inter.coordinates={[0.00 0.00 0.00]`.
- Lines 27: computes `inter.relaxation` using `inter.relaxation={'damp'}`.
- Lines 28: computes `inter.rlx_keep` using `inter.rlx_keep='diagonal'`.
- Lines 29: computes `inter.equilibrium` using `inter.equilibrium='zero'`.
- Lines 30: computes `inter.damp_rate` using `inter.damp_rate=300`.
- Lines 33: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 34: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 37: computes `sys.disable` using `sys.disable={'krylov','trajlevel'}`.
- Lines 40: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 44: computes `theta` using `theta=atan(sqrt(2))`.
- Lines 47: computes `parameters.max_rank` using `parameters.max_rank=7`.
- Lines 48: computes `parameters.axis` using `parameters.axis=[sqrt(2/3) 0 sqrt(1/3)]`.
- Lines 49: computes `parameters.rate` using `parameters.rate=-19840`.
- Lines 50: computes `parameters.grid` using `parameters.grid='rep_2ang_200pts_oct'`.

## Implementation structure

- Cross-polarization experiment between protons and 14N overtone
- transition in glycine under MAS. Glycine quadrupolar tensor da-
- ta comes from the paper by O'Dell and Ratcliffe:
- Magnetisation accumulation profile as a function of contact ti-
- me with a rough powder grid.
- Calculation time: minutes
- System specification
- Relaxation theory
- Basis set
- Algorithmic options
- Spinach housekeeping
- Magic angle

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `eeqq2nqi()`, `create()`, `basis()`, `atan()`, `state()`, `operator()`, `kfigure()`, `scale_figure()`, `subplot()`, `singlerot()`, `plot_1d()`, `set()`, `ktitle()`, `num2str()`, `kxlabel()`.
