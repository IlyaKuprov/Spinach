# examples/esr_liq_pulsed/endor_phenyl.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/esr_liq_pulsed/endor_phenyl.m`
- Signature: `endor_phenyl()`
- Total lines: 69

## Purpose

Mims ENDOR on a phenyl radical in liquid state. The g-factor and the isotropic proton hyperfine couplings are specified explicitly from the experimental data of Kasai, Hedaya, and Whipple (J. Am. Chem. Soc. 1969, 91, 4364): a(ortho)=17.4 G, a(meta)=5.9 G, a(para)=1.9 G, and an isotropic g-factor of 2.0024. The two ortho and the two meta protons are magnetically equivalent, so the full symmetry treatment uses an S2 x 

## Physical / mathematical content

- Liquid-state ESR examples. The dominant physics is electron Zeeman interaction, hyperfine coupling, relaxation broadening, and pulse-acquire or ENDOR-type detection in fast tumbling systems.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 16-17: Magnet field; implemented by `sys.magnet=0.33`.
- Lines 19-20: Electron and the five ring protons (two ortho, two meta, one para); implemented by `sys.isotopes={'E','1H','1H','1H','1H','1H'}`.
- Lines 22-23: Phenyl radical isotropic g-factor; implemented by `g_phenyl=2.0024`.
- Lines 26-28: Isotropic proton hyperfine couplings, converted from milliTesla using the phenyl radical g-factor; implemented by `inter.coupling.scalar=cell(6,6)`.
- Lines 35-36: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 39-40: Symmetry: equivalent ortho pair and equivalent meta pair; implemented by `bas.sym_group={'S2','S2'}`.
- Lines 43-44: Sequence parameters; implemented by `parameters.offset=0`.
- Lines 52-53: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 56-57: Simulation; implemented by `fid=liquid(spin_system,@endor_mims,parameters,'esr')`.
- Lines 59-60: Crude apodisation; implemented by `fid=apodisation(spin_system,fid-mean(fid),{{'kaiser',6}})`.
- Lines 62-63: Fourier transform; implemented by `spectrum=fftshift(fft(fid,parameters.zerofill))`.
- Lines 65-66: Plotting; implemented by `kfigure(); plot_1d(spin_system,abs(spectrum),parameters)`.

### Key state/data transformations

- Lines 17: computes `sys.magnet` using `sys.magnet=0.33`.
- Lines 20: computes `sys.isotopes` using `sys.isotopes={'E','1H','1H','1H','1H','1H'}`.
- Lines 23: computes `g_phenyl` using `g_phenyl=2.0024`.
- Lines 24: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={g_phenyl 0 0 0 0 0}`.
- Lines 28: computes `inter.coupling.scalar` using `inter.coupling.scalar=cell(6,6)`.
- Lines 29: computes `inter.coupling.scalar{2,1}` using `inter.coupling.scalar{2,1}=mt2hz(1.74,g_phenyl)`.
- Lines 30: computes `inter.coupling.scalar{3,1}` using `inter.coupling.scalar{3,1}=mt2hz(1.74,g_phenyl)`.
- Lines 31: computes `inter.coupling.scalar{4,1}` using `inter.coupling.scalar{4,1}=mt2hz(0.59,g_phenyl)`.
- Lines 32: computes `inter.coupling.scalar{5,1}` using `inter.coupling.scalar{5,1}=mt2hz(0.59,g_phenyl)`.
- Lines 33: computes `inter.coupling.scalar{6,1}` using `inter.coupling.scalar{6,1}=mt2hz(0.19,g_phenyl)`.
- Lines 36: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 37: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 40: computes `bas.sym_group` using `bas.sym_group={'S2','S2'}`.
- Lines 41: computes `bas.sym_spins` using `bas.sym_spins={[2 3],[4 5]}`.
- Lines 44: computes `parameters.offset` using `parameters.offset=0`.
- Lines 45: computes `parameters.npoints` using `parameters.npoints=512`.
- Lines 46: computes `parameters.sweep` using `parameters.sweep=3e8`.
- Lines 47: computes `parameters.tau` using `parameters.tau=100e-9`.

## Implementation structure

- Mims ENDOR on a phenyl radical in liquid state. The g-factor and the
- isotropic proton hyperfine couplings are specified explicitly from the
- experimental data of Kasai, Hedaya, and Whipple (J. Am. Chem. Soc.
- 1969, 91, 4364): a(ortho)=17.4 G, a(meta)=5.9 G, a(para)=1.9 G, and an
- isotropic g-factor of 2.0024. The two ortho and the two meta protons
- are magnetically equivalent, so the full symmetry treatment uses an
- S2 x S2 group direct product.
- Calculation time: seconds
- Magnet field
- Electron and the five ring protons (two ortho, two meta, one para)
- Phenyl radical isotropic g-factor
- Isotropic proton hyperfine couplings, converted from milliTesla

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `mt2hz()`, `create()`, `basis()`, `liquid()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
