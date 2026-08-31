# examples/nmr_solids/static_powder_nqi_a.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_solids/static_powder_nqi_a.m`
- Signature: `static_powder_nqi_a()`
- Total lines: 54

## Purpose

Static quadrupolar 14N powder pattern of L-valyl-L-alanine using very large numerical orientation grid, set to reproduce Figure 5 from the paper by O'Dell and Ratcliffe: Calculation time: minutes

## Physical / mathematical content

- Solid-state NMR examples. The key physics is anisotropic spin interactions under static or magic-angle-spinning conditions: chemical-shift anisotropy, dipolar coupling, quadrupolar coupling, cross-polarisation, and orientation averaging using Floquet, Fokker-Planck, or direct powder quadrature formalisms.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.
- Quadrupolar physics is relevant: nuclei with spin > 1/2 interact with the electric field gradient tensor, introducing second-rank anisotropy, asymmetry, and overtone or MQ phenomena.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 13-14: System specification; implemented by `sys.magnet=21.1`.
- Lines 19-20: Basis set; implemented by `bas.formalism='zeeman-hilb'`.
- Lines 23-24: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 27-28: Sequence parameters; implemented by `parameters.spins={'14N'}`.
- Lines 41-42: Simulation; implemented by `fid=powder(spin_system,@acquire,parameters,'nmr')`.
- Lines 44-45: Apodisation; implemented by `fid=apodisation(spin_system,fid,{{'exp',6}})`.
- Lines 47-48: Fourier transform; implemented by `spectrum=fftshift(fft(fid,parameters.zerofill))`.
- Lines 50-51: Plotting; implemented by `kfigure(); plot_1d(spin_system,real(spectrum),parameters)`.

### Key state/data transformations

- Lines 14: computes `sys.magnet` using `sys.magnet=21.1`.
- Lines 15: computes `sys.isotopes` using `sys.isotopes={'14N','14N'}`.
- Lines 16: computes `inter.coupling.matrix{1,1}` using `inter.coupling.matrix{1,1}=eeqq2nqi(1.24e6,0.22,1,[0 0 0])`.
- Lines 17: computes `inter.coupling.matrix{2,2}` using `inter.coupling.matrix{2,2}=eeqq2nqi(3.06e6,0.40,1,[0 0 0])`.
- Lines 20: computes `bas.formalism` using `bas.formalism='zeeman-hilb'`.
- Lines 21: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 24: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 28: computes `parameters.spins` using `parameters.spins={'14N'}`.
- Lines 29: computes `parameters.decouple` using `parameters.decouple={}`.
- Lines 30: computes `parameters.offset` using `parameters.offset=0`.
- Lines 31: computes `parameters.sweep` using `parameters.sweep=6e6`.
- Lines 32: computes `parameters.npoints` using `parameters.npoints=512`.
- Lines 33: computes `parameters.zerofill` using `parameters.zerofill=2048`.
- Lines 34: computes `parameters.axis_units` using `parameters.axis_units='MHz'`.
- Lines 35: computes `parameters.invert_axis` using `parameters.invert_axis=1`.
- Lines 36: computes `parameters.grid` using `parameters.grid='icos_2ang_163842pts'`.
- Lines 37: computes `parameters.rho0` using `parameters.rho0=state(spin_system,'L+','14N')`.
- Lines 38: computes `parameters.coil` using `parameters.coil=state(spin_system,'L+','14N')`.

## Implementation structure

- Static quadrupolar 14N powder pattern of L-valyl-L-alanine using
- very large numerical orientation grid, set to reproduce Figure 5
- from the paper by O'Dell and Ratcliffe:
- Calculation time: minutes
- System specification
- Basis set
- Spinach housekeeping
- Sequence parameters
- Simulation
- Apodisation
- Fourier transform
- Plotting

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `eeqq2nqi()`, `create()`, `basis()`, `state()`, `powder()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
