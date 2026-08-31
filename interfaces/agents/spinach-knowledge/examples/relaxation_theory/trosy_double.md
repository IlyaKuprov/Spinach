# examples/relaxation_theory/trosy_double.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/relaxation_theory/trosy_double.m`
- Signature: `trosy_double()`
- Total lines: 85

## Purpose

Hari Arthanari's Double TROSY effect. Calculation time: seconds.

## Physical / mathematical content

- Relaxation-theory examples. The mathematical backbone is Bloch-Redfield-Wangsness or stochastic Liouville theory, spectral densities, cross-correlation terms, motional models, and extraction of longitudinal/transverse decay behaviour from superoperators.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.
- Chemical-shift anisotropy is present: shielding is treated as a second-rank tensor whose orientation relative to the field or rotor axis modulates line shapes and transfer dynamics.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 9-10: Magnet field; implemented by `sys.magnet=14.1`.
- Lines 12-14: Read 3-fluorotyrosine DFT calculation; implemented by `[~,inter_dft]=g2spinach(gparse('../standard_systems/4_fluoro_phe.out'), {{'H','1H'},{'C','13C'},{'F','19F'}},[32.07 186.38 192.97])`.
- Lines 16-17: Extract coordinates and CSAs; implemented by `sys.isotopes={'1H','1H','19F','13C'}`.
- Lines 29-30: Extract J-couplings; implemented by `inter.coupling.scalar=inter_dft.coupling.scalar([10 20 19 8],[10 20 19 8])`.
- Lines 32-33: Relaxation theory; implemented by `inter.relaxation={'redfield'}`.
- Lines 38-39: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 42-43: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 46-47: Sequence parameters -13C; implemented by `parameters.spins={'13C'}`.
- Lines 58-59: Simulation; implemented by `fid=liquid(spin_system,@acquire,parameters,'nmr')`.
- Lines 61-62: Apodisation; implemented by `fid=apodisation(spin_system,fid,{{'gauss',6}})`.
- Lines 64-65: Fourier transform; implemented by `spectrum=fftshift(fft(fid,parameters.zerofill))`.
- Lines 67-68: Plotting, full spectrum; implemented by `kfigure(); subplot(2,1,1)`.
- Lines 72-73: Re-run with no DD to protons; implemented by `inter.coordinates{1}=[]`.

### Key state/data transformations

- Lines 10: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 13-14: computes `[~,inter_dft]` using `[~,inter_dft]=g2spinach(gparse('../standard_systems/4_fluoro_phe.out'), {{'H','1H'},{'C','13C'},{'F','19F'}},[32.07 186.38 192.97])`.
- Lines 17: computes `sys.isotopes` using `sys.isotopes={'1H','1H','19F','13C'}`.
- Lines 18: computes `inter.zeeman.matrix` using `inter.zeeman.matrix=cell(1,4)`.
- Lines 19: computes `inter.zeeman.matrix{1}` using `inter.zeeman.matrix{1}=inter_dft.zeeman.matrix{10}`.
- Lines 20: computes `inter.zeeman.matrix{2}` using `inter.zeeman.matrix{2}=inter_dft.zeeman.matrix{20}`.
- Lines 21: computes `inter.zeeman.matrix{3}` using `inter.zeeman.matrix{3}=inter_dft.zeeman.matrix{19}`.
- Lines 22: computes `inter.zeeman.matrix{4}` using `inter.zeeman.matrix{4}=inter_dft.zeeman.matrix{8}`.
- Lines 23: computes `inter.coordinates` using `inter.coordinates=cell(4,1)`.
- Lines 24: computes `inter.coordinates{1}` using `inter.coordinates{1}=inter_dft.coordinates{10}`.
- Lines 25: computes `inter.coordinates{2}` using `inter.coordinates{2}=inter_dft.coordinates{20}`.
- Lines 26: computes `inter.coordinates{3}` using `inter.coordinates{3}=inter_dft.coordinates{19}`.
- Lines 27: computes `inter.coordinates{4}` using `inter.coordinates{4}=inter_dft.coordinates{8}`.
- Lines 30: computes `inter.coupling.scalar` using `inter.coupling.scalar=inter_dft.coupling.scalar([10 20 19 8],[10 20 19 8])`.
- Lines 33: computes `inter.relaxation` using `inter.relaxation={'redfield'}`.
- Lines 34: computes `inter.rlx_keep` using `inter.rlx_keep='secular'`.
- Lines 35: computes `inter.equilibrium` using `inter.equilibrium='zero'`.
- Lines 36: computes `inter.tau_c` using `inter.tau_c={20e-9}`.

## Implementation structure

- Hari Arthanari's Double TROSY effect.
- Calculation time: seconds.
- Magnet field
- Read 3-fluorotyrosine DFT calculation
- Extract coordinates and CSAs
- Extract J-couplings
- Relaxation theory
- Basis set
- Spinach housekeeping
- Sequence parameters -13C
- Simulation
- Apodisation

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `g2spinach()`, `gparse()`, `create()`, `basis()`, `state()`, `liquid()`, `apodisation()`, `fftshift()`, `kfigure()`, `subplot()`, `plot_1d()`, `kxlabel()`.
