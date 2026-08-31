# examples/nmr_solids/case_studies/akbey_2h_13c_mas/fig2_three_site.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_solids/case_studies/akbey_2h_13c_mas/fig2_three_site.m`
- Signature: `fig2_three_site()`
- Total lines: 72

## Purpose

Three-site position exchange for a deuterium nucleus. The sites differ in the chemical shift and the orientation of the quadru- polar tensor. Set to reproduce Figure 2 in: Calculation time: seconds.

## Physical / mathematical content

- Solid-state NMR examples. The key physics is anisotropic spin interactions under static or magic-angle-spinning conditions: chemical-shift anisotropy, dipolar coupling, quadrupolar coupling, cross-polarisation, and orientation averaging using Floquet, Fokker-Planck, or direct powder quadrature formalisms.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.
- Quadrupolar physics is relevant: nuclei with spin > 1/2 interact with the electric field gradient tensor, introducing second-rank anisotropy, asymmetry, and overtone or MQ phenomena.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 14-15: Magnet field; implemented by `sys.magnet=9.4`.
- Lines 17-18: Spin system; implemented by `sys.isotopes={'2H','2H', '2H'}`.
- Lines 20-21: Quadrupolar interactions; implemented by `[Q1,Q2, Q3]=weblab2nqi(0.16e6,0,1,0,acos(1/3),2*pi/3)`.
- Lines 27-28: Kinetics; implemented by `inter.chem.rates=1e4*[-2 1 1`.
- Lines 33-34: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 37-38: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 41-42: Sequence parameters; implemented by `parameters.spins={'2H'}`.
- Lines 54-55: MAS parameters; implemented by `parameters.rate=8500`.
- Lines 59-60: Acquisition; implemented by `fid=singlerot(spin_system,@acquire,parameters,'nmr')`.
- Lines 62-63: Apodisation; implemented by `fid=apodisation(spin_system,fid,{{'none'}})`.
- Lines 65-66: Fourier transform; implemented by `spectrum=fftshift(fft(fid,parameters.zerofill))`.
- Lines 68-69: Plotting; implemented by `kfigure(); plot_1d(spin_system,real(spectrum),parameters)`.

### Key state/data transformations

- Lines 15: computes `sys.magnet` using `sys.magnet=9.4`.
- Lines 18: computes `sys.isotopes` using `sys.isotopes={'2H','2H', '2H'}`.
- Lines 21: computes `[Q1,Q2, Q3]` using `[Q1,Q2, Q3]=weblab2nqi(0.16e6,0,1,0,acos(1/3),2*pi/3)`.
- Lines 22: computes `inter.coupling.matrix{1,1}` using `inter.coupling.matrix{1,1}=Q1`.
- Lines 23: computes `inter.coupling.matrix{2,2}` using `inter.coupling.matrix{2,2}=Q2`.
- Lines 24: computes `inter.coupling.matrix{3,3}` using `inter.coupling.matrix{3,3}=Q3`.
- Lines 25: computes `inter.chem.parts` using `inter.chem.parts={1,2,3}`.
- Lines 28: computes `inter.chem.rates` using `inter.chem.rates=1e4*[-2 1 1`.
- Lines 31: computes `inter.chem.concs` using `inter.chem.concs=[1.0 1.0 1.0]`.
- Lines 34: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 35: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 38: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 42: computes `parameters.spins` using `parameters.spins={'2H'}`.
- Lines 43: computes `parameters.rho0` using `parameters.rho0=state(spin_system,'L+','2H','chem')`.
- Lines 44: computes `parameters.coil` using `parameters.coil=state(spin_system,'L+','2H')`.
- Lines 45: computes `parameters.decouple` using `parameters.decouple={}`.
- Lines 46: computes `parameters.offset` using `parameters.offset=0`.
- Lines 47: computes `parameters.sweep` using `parameters.sweep=0.4e6`.

## Implementation structure

- Three-site position exchange for a deuterium nucleus. The sites
- differ in the chemical shift and the orientation of the quadru-
- polar tensor. Set to reproduce Figure 2 in:
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
