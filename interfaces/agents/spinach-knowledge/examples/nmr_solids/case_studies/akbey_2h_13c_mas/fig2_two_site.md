# examples/nmr_solids/case_studies/akbey_2h_13c_mas/fig2_two_site.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_solids/case_studies/akbey_2h_13c_mas/fig2_two_site.m`
- Signature: `fig2_two_site()`
- Total lines: 70

## Purpose

Two-site position exchange for a deuterium nucleus. The sites differ in the chemical shift and the orientation of the quad- rupolar tensor. Set to reproduce Figure 2 in: Calculation time: seconds.

## Physical / mathematical content

- Solid-state NMR examples. The key physics is anisotropic spin interactions under static or magic-angle-spinning conditions: chemical-shift anisotropy, dipolar coupling, quadrupolar coupling, cross-polarisation, and orientation averaging using Floquet, Fokker-Planck, or direct powder quadrature formalisms.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.
- Quadrupolar physics is relevant: nuclei with spin > 1/2 interact with the electric field gradient tensor, introducing second-rank anisotropy, asymmetry, and overtone or MQ phenomena.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 14-15: Magnet field; implemented by `sys.magnet=9.4`.
- Lines 17-18: Spin system; implemented by `sys.isotopes={'2H','2H'}`.
- Lines 20-21: Quadrupolar interactions; implemented by `[Q1,Q2]=weblab2nqi(0.16e6,0,1,0,acos(1/3),2*pi/3)`.
- Lines 25-26: Kinetics; implemented by `inter.chem.parts={1,2}`.
- Lines 31-32: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 35-36: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 39-40: Sequence parameters; implemented by `parameters.spins={'2H'}`.
- Lines 52-53: MAS parameters; implemented by `parameters.rate=8500`.
- Lines 57-58: Acquisition; implemented by `fid=singlerot(spin_system,@acquire,parameters,'nmr')`.
- Lines 60-61: Apodisation; implemented by `fid=apodisation(spin_system,fid,{{'none'}})`.
- Lines 63-64: Fourier transform; implemented by `spectrum=fftshift(fft(fid,parameters.zerofill))`.
- Lines 66-67: Plotting; implemented by `kfigure(); plot_1d(spin_system,real(spectrum),parameters)`.

### Key state/data transformations

- Lines 15: computes `sys.magnet` using `sys.magnet=9.4`.
- Lines 18: computes `sys.isotopes` using `sys.isotopes={'2H','2H'}`.
- Lines 21: computes `[Q1,Q2]` using `[Q1,Q2]=weblab2nqi(0.16e6,0,1,0,acos(1/3),2*pi/3)`.
- Lines 22: computes `inter.coupling.matrix{1,1}` using `inter.coupling.matrix{1,1}=Q1`.
- Lines 23: computes `inter.coupling.matrix{2,2}` using `inter.coupling.matrix{2,2}=Q2`.
- Lines 26: computes `inter.chem.parts` using `inter.chem.parts={1,2}`.
- Lines 27: computes `inter.chem.rates` using `inter.chem.rates=1e4*[-1 1`.
- Lines 29: computes `inter.chem.concs` using `inter.chem.concs=[1.0 1.0]`.
- Lines 32: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 33: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 36: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 40: computes `parameters.spins` using `parameters.spins={'2H'}`.
- Lines 41: computes `parameters.rho0` using `parameters.rho0=state(spin_system,'L+','2H','chem')`.
- Lines 42: computes `parameters.coil` using `parameters.coil=state(spin_system,'L+','2H')`.
- Lines 43: computes `parameters.decouple` using `parameters.decouple={}`.
- Lines 44: computes `parameters.offset` using `parameters.offset=0`.
- Lines 45: computes `parameters.sweep` using `parameters.sweep=0.4e6`.
- Lines 46: computes `parameters.npoints` using `parameters.npoints=1024`.

## Implementation structure

- Two-site position exchange for a deuterium nucleus. The sites
- differ in the chemical shift and the orientation of the quad-
- rupolar tensor. Set to reproduce Figure 2 in:
- Calculation time: seconds.
- Magnet field
- Spin system
- Quadrupolar interactions
- Kinetics
- Basis set
- Spinach housekeeping
- Sequence parameters
- MAS parameters

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `weblab2nqi()`, `acos()`, `create()`, `basis()`, `state()`, `singlerot()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
