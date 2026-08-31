# examples/nmr_overtone/cpmas_glycine_match_1.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_overtone/cpmas_glycine_match_1.m`
- Signature: `cpmas_glycine_match_1()`
- Total lines: 92

## Purpose

Cross-polarization experiment between protons and 14N overtone transition in glycine under MAS. Glycine quadrupolar tensor da- ta comes from the paper by O'Dell and Ratcliffe: Hartmann-Hahn condition profile with a rough powder grid, as a function of 1H RF power. Calculation time: minutes

## Physical / mathematical content

- Overtone NMR examples. The important regime is excitation or detection of formally forbidden high-order transitions in quadrupolar nuclei, usually aided by MAS or Fokker-Planck treatments of periodic motion.
- Quadrupolar physics is relevant: nuclei with spin > 1/2 interact with the electric field gradient tensor, introducing second-rank anisotropy, asymmetry, and overtone or MQ phenomena.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 18-19: System specification; implemented by `sys.magnet=14.10220742; sys.isotopes={'14N','1H'}`.
- Lines 26-27: Relaxation theory; implemented by `inter.relaxation={'damp'}`.
- Lines 32-33: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 36-37: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 40-41: Magic angle; implemented by `theta=atan(sqrt(2))`.
- Lines 43-44: Spectrum setup; implemented by `parameters.max_rank=7`.
- Lines 65-66: Proton RF power array; implemented by `rf_powers=linspace(25e3,39e3,15)`.
- Lines 68-69: Start a new figure; implemented by `kfigure(); scale_figure([3.0 1.0])`.
- Lines 71-72: Proton RF power scan; implemented by `for n=1:15`.
- Lines 74-75: Subplot selection; implemented by `subplot(1,15,n)`.
- Lines 77-78: Power setting; implemented by `parameters.rf_pwr=2*pi*[55e3 rf_powers(n)]/sin(theta)`.
- Lines 80-81: Simulation; implemented by `spectrum=singlerot(spin_system,@overtone_cp,parameters,'qnmr')`.
- Lines 83-84: Plotting; implemented by `plot_1d(spin_system,real(spectrum),parameters)`.

### Control flow inferred from the code

- Line 72: `for` loop over `n=1:15`.

### Key state/data transformations

- Lines 19: computes `sys.magnet` using `sys.magnet=14.10220742; sys.isotopes={'14N','1H'}`.
- Lines 20: computes `inter.coupling.matrix{1,1}` using `inter.coupling.matrix{1,1}=eeqq2nqi(1.18e6,0.53,1,[0 0 0])`.
- Lines 21: computes `inter.coupling.matrix{2,2}` using `inter.coupling.matrix{2,2}=[]`.
- Lines 22: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={32.4 0.0}`.
- Lines 23: computes `inter.coordinates` using `inter.coordinates={[0.00 0.00 0.00]`.
- Lines 27: computes `inter.relaxation` using `inter.relaxation={'damp'}`.
- Lines 28: computes `inter.rlx_keep` using `inter.rlx_keep='diagonal'`.
- Lines 29: computes `inter.equilibrium` using `inter.equilibrium='zero'`.
- Lines 30: computes `inter.damp_rate` using `inter.damp_rate=300`.
- Lines 33: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 34: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 37: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 41: computes `theta` using `theta=atan(sqrt(2))`.
- Lines 44: computes `parameters.max_rank` using `parameters.max_rank=7`.
- Lines 45: computes `parameters.axis` using `parameters.axis=[sqrt(2/3) 0 sqrt(1/3)]`.
- Lines 46: computes `parameters.rate` using `parameters.rate=-19840`.
- Lines 47: computes `parameters.grid` using `parameters.grid='rep_2ang_200pts_oct'`.
- Lines 48: computes `parameters.sweep` using `parameters.sweep=[4.4e4 5.2e4]`.

## Implementation structure

- Cross-polarization experiment between protons and 14N overtone
- transition in glycine under MAS. Glycine quadrupolar tensor da-
- ta comes from the paper by O'Dell and Ratcliffe:
- Hartmann-Hahn condition profile with a rough powder grid, as a
- function of 1H RF power.
- Calculation time: minutes
- System specification
- Relaxation theory
- Basis set
- Spinach housekeeping
- Magic angle
- Spectrum setup

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `eeqq2nqi()`, `create()`, `basis()`, `atan()`, `state()`, `operator()`, `kfigure()`, `scale_figure()`, `subplot()`, `rf_powers()`, `singlerot()`, `plot_1d()`, `set()`, `ktitle()`, `num2str()`, `kxlabel()`.
