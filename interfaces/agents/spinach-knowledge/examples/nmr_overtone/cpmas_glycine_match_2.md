# examples/nmr_overtone/cpmas_glycine_match_2.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_overtone/cpmas_glycine_match_2.m`
- Signature: `cpmas_glycine_match_2()`
- Total lines: 110

## Purpose

Cross-polarization experiment between protons and 14N overtone transition in glycine under MAS. Glycine quadrupolar tensor da- ta comes from the paper by O'Dell and Ratcliffe: Hartmann-Hahn condition profile with a rough powder grid, as a function of spinning rate and 1H RF power. Calculation time: hours

## Physical / mathematical content

- Overtone NMR examples. The important regime is excitation or detection of formally forbidden high-order transitions in quadrupolar nuclei, usually aided by MAS or Fokker-Planck treatments of periodic motion.
- Quadrupolar physics is relevant: nuclei with spin > 1/2 interact with the electric field gradient tensor, introducing second-rank anisotropy, asymmetry, and overtone or MQ phenomena.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 18-19: System specification; implemented by `sys.magnet=14.10220742; sys.isotopes={'14N','1H'}`.
- Lines 26-27: Relaxation theory; implemented by `inter.relaxation={'damp'}`.
- Lines 32-33: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 36-37: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 40-41: Magic angle; implemented by `theta=atan(sqrt(2))`.
- Lines 43-44: Spectrum setup; implemented by `parameters.max_rank=5`.
- Lines 62-63: Proton RF power array; implemented by `rf_powers=linspace(10e3,200e3,50)`.
- Lines 65-66: Spinning rate array; implemented by `spin_rates=linspace(20e3,90e3,50)`.
- Lines 68-69: Matrix preallocation; implemented by `intensities=zeros(numel(rf_powers),numel(spin_rates))`.
- Lines 71-72: Loop over powers and spinning rates; implemented by `parfor nk=1:numel(rf_powers)*numel(spin_rates)`.
- Lines 74-75: Get the indices; implemented by `[n,k]=ind2sub([numel(rf_powers) numel(spin_rates)],nk)`.
- Lines 77-78: Get a copy of parameters; implemented by `localpar=parameters`.
- Lines 80-81: Power setting; implemented by `localpar.rf_pwr=2*pi*[55e3 rf_powers(n)]/sin(theta)`.
- Lines 83-84: Spinning rate; implemented by `localpar.rate=-spin_rates(k)`.
- Lines 86-87: RF frequency on nitrogen OT; implemented by `localpar.rf_frq=8e3-2*localpar.rate`.
- Lines 89-91: Sweep width; implemented by `localpar.sweep=[localpar.rf_frq-4e3 localpar.rf_frq+4e3]`.
- Lines 93-94: Simulation; implemented by `spectrum=singlerot(spin_system,@overtone_cp,localpar,'qnmr')`.
- Lines 96-97: Intensity; implemented by `intensities(nk)=sum(real(spectrum))`.

### Control flow inferred from the code

- Line 72: `parfor` loop over `nk=1:numel(rf_powers)*numel(spin_rates)`.

### Key state/data transformations

- Lines 19: computes `sys.magnet` using `sys.magnet=14.10220742; sys.isotopes={'14N','1H'}`.
- Lines 20: computes `inter.coupling.matrix{1,1}` using `inter.coupling.matrix{1,1}=eeqq2nqi(1.18e6,0.53,1,[0 0 0])`.
- Lines 21: computes `inter.coupling.matrix{2,2}` using `inter.coupling.matrix{2,2}=[]`.
- Lines 22: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={32.4 0.0}`.
- Lines 23: computes `inter.coordinates` using `inter.coordinates={[0.00 0.00 0.00]`.
- Lines 27: computes `inter.relaxation` using `inter.relaxation={'damp'}`.
- Lines 28: computes `inter.rlx_keep` using `inter.rlx_keep='diagonal'`.
- Lines 29: computes `inter.equilibrium` using `inter.equilibrium='zero'`.
- Lines 30: computes `inter.damp_rate` using `inter.damp_rate=1000`.
- Lines 33: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 34: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 37: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 41: computes `theta` using `theta=atan(sqrt(2))`.
- Lines 44: computes `parameters.max_rank` using `parameters.max_rank=5`.
- Lines 45: computes `parameters.axis` using `parameters.axis=[sqrt(2/3) 0 sqrt(1/3)]`.
- Lines 46: computes `parameters.grid` using `parameters.grid='rep_2ang_200pts_oct'`.
- Lines 47: computes `parameters.npoints` using `parameters.npoints=256`.
- Lines 48: computes `parameters.zerofill` using `parameters.zerofill=256`.

## Implementation structure

- Cross-polarization experiment between protons and 14N overtone
- transition in glycine under MAS. Glycine quadrupolar tensor da-
- ta comes from the paper by O'Dell and Ratcliffe:
- Hartmann-Hahn condition profile with a rough powder grid, as
- a function of spinning rate and 1H RF power.
- Calculation time: hours
- System specification
- Relaxation theory
- Basis set
- Spinach housekeeping
- Magic angle
- Spectrum setup

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `eeqq2nqi()`, `create()`, `basis()`, `atan()`, `state()`, `operator()`, `ind2sub()`, `rf_powers()`, `spin_rates()`, `singlerot()`, `intensities()`, `kfigure()`, `kxlabel()`, `kylabel()`, `set()`.
