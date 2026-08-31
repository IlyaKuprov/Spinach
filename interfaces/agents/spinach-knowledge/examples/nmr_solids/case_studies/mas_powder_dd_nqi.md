# examples/nmr_solids/case_studies/mas_powder_dd_nqi.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_solids/case_studies/mas_powder_dd_nqi.m`
- Signature: `mas_powder_dd_nqi()`
- Total lines: 67

## Purpose

Powder magic angle spinning spectrum of a pair of dipole-coupled quadrupolar nuclei; this is apparently something that other simu- lation packages cannot do. Parameters from Jeongjae Lee. Calculation time: seconds

## Physical / mathematical content

- Solid-state NMR examples. The key physics is anisotropic spin interactions under static or magic-angle-spinning conditions: chemical-shift anisotropy, dipolar coupling, quadrupolar coupling, cross-polarisation, and orientation averaging using Floquet, Fokker-Planck, or direct powder quadrature formalisms.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.
- Quadrupolar physics is relevant: nuclei with spin > 1/2 interact with the electric field gradient tensor, introducing second-rank anisotropy, asymmetry, and overtone or MQ phenomena.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 11-12: System specification; implemented by `sys.magnet=9.4`.
- Lines 15-16: Interactions; implemented by `inter.coordinates={[2.602 8.750 3.651]`.
- Lines 25-26: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 30-34: Enable GPU sys.enable={'gpu'};; implemented by `spin_system=create(sys,inter)`.
- Lines 33-34: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 37-38: Experiment setup; implemented by `parameters.rate=100000`.
- Lines 54-55: Simulation; implemented by `fid=singlerot(spin_system,@acquire,parameters,'nmr')`.
- Lines 57-58: Apodisation; implemented by `fid=apodisation(spin_system,fid,{{'exp',6}})`.
- Lines 60-61: Fourier transform; implemented by `spectrum=fftshift(fft(fid,parameters.zerofill))`.
- Lines 63-64: Plotting; implemented by `kfigure(); plot_1d(spin_system,real(spectrum),parameters)`.

### Key state/data transformations

- Lines 12: computes `sys.magnet` using `sys.magnet=9.4`.
- Lines 13: computes `sys.isotopes` using `sys.isotopes={'23Na','17O'}`.
- Lines 16: computes `inter.coordinates` using `inter.coordinates={[2.602 8.750 3.651]`.
- Lines 18: computes `inter.coupling.matrix{1,1}` using `inter.coupling.matrix{1,1}=castep2nqi([-0.0497 0.0520 -0.0019`.
- Lines 21: computes `inter.coupling.matrix{2,2}` using `inter.coupling.matrix{2,2}=castep2nqi([ 0.1580 0.0340 -0.5562`.
- Lines 26: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 27: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 28: computes `bas.projections` using `bas.projections=+1`.
- Lines 34: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 38: computes `parameters.rate` using `parameters.rate=100000`.
- Lines 39: computes `parameters.axis` using `parameters.axis=[sqrt(2/3) 0 sqrt(1/3)]`.
- Lines 40: computes `parameters.max_rank` using `parameters.max_rank=30`.
- Lines 41: computes `parameters.grid` using `parameters.grid='rep_2ang_200pts_sph'`.
- Lines 42: computes `parameters.sweep` using `parameters.sweep=5e6`.
- Lines 43: computes `parameters.npoints` using `parameters.npoints=1024`.
- Lines 44: computes `parameters.zerofill` using `parameters.zerofill=4096`.
- Lines 45: computes `parameters.offset` using `parameters.offset=0`.
- Lines 46: computes `parameters.spins` using `parameters.spins={'17O'}`.

## Implementation structure

- Powder magic angle spinning spectrum of a pair of dipole-coupled
- quadrupolar nuclei; this is apparently something that other simu-
- lation packages cannot do. Parameters from Jeongjae Lee.
- Calculation time: seconds
- System specification
- Interactions
- Basis set
- Enable GPU
- sys.enable={'gpu'};
- Spinach housekeeping
- Experiment setup
- Simulation

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `castep2nqi()`, `create()`, `basis()`, `state()`, `singlerot()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
