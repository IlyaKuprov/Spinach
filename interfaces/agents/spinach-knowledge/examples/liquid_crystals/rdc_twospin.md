# examples/liquid_crystals/rdc_twospin.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/liquid_crystals/rdc_twospin.m`
- Signature: `rdc_twospin()`
- Total lines: 64

## Purpose

CLIP-HSQC spectrum of a C-H system in a liquid crystal with a user-specified order matrix. Calculation time: seconds.

## Physical / mathematical content

- Liquid-crystal examples. These scripts exploit partial ordering and Saupe-tensor physics, so anisotropic couplings survive orientational averaging and generate residual dipolar couplings or anisotropic transfer behaviour.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 10-11: Spin system parameters; implemented by `sys.magnet=5.9`.
- Lines 20-21: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 24-25: Algorithmic options; implemented by `sys.disable={'colorbar'}`.
- Lines 27-28: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 31-32: Sequence parameters; implemented by `parameters.sweep=[3000 1000]`.
- Lines 41-42: Simulation; implemented by `fid=liquid(spin_system,@clip_hsqc,parameters,'nmr')`.
- Lines 44-45: Apodisation; implemented by `fid.pos=apodisation(spin_system,fid.pos,{{'sqcos'},{'sqcos'}})`.
- Lines 48-49: F2 Fourier transform; implemented by `f1_pos=fftshift(fft(fid.pos,parameters.zerofill(2),1),1)`.
- Lines 52-53: Form States signal; implemented by `fid=f1_pos+conj(f1_neg)`.
- Lines 55-56: F1 Fourier transform; implemented by `spectrum=fftshift(fft(fid,parameters.zerofill(1),2),2)`.
- Lines 58-59: Plotting; implemented by `kfigure(); scale_figure([1.5 2.0])`.

### Key state/data transformations

- Lines 11: computes `sys.magnet` using `sys.magnet=5.9`.
- Lines 12: computes `sys.isotopes` using `sys.isotopes={'1H','13C'}`.
- Lines 13: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={5.0 65.0}`.
- Lines 14: computes `inter.coupling.scalar` using `inter.coupling.scalar=cell(2)`.
- Lines 15: computes `inter.coupling.scalar{1,2}` using `inter.coupling.scalar{1,2}=140`.
- Lines 16: computes `inter.coordinates` using `inter.coordinates={[0.0 0.0 0.0]`.
- Lines 18: computes `inter.order_matrix` using `inter.order_matrix={diag([1e-3 2e-3 -3e-3])}`.
- Lines 21: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 22: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 25: computes `sys.disable` using `sys.disable={'colorbar'}`.
- Lines 28: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 32: computes `parameters.sweep` using `parameters.sweep=[3000 1000]`.
- Lines 33: computes `parameters.offset` using `parameters.offset=[4250 1200]`.
- Lines 34: computes `parameters.npoints` using `parameters.npoints=[128 128]`.
- Lines 35: computes `parameters.zerofill` using `parameters.zerofill=[512 512]`.
- Lines 36: computes `parameters.spins` using `parameters.spins={'13C','1H'}`.
- Lines 37: computes `parameters.axis_units` using `parameters.axis_units='ppm'`.
- Lines 38: computes `parameters.needs` using `parameters.needs={'rdc'}`.

## Implementation structure

- CLIP-HSQC spectrum of a C-H system in a liquid crystal
- with a user-specified order matrix.
- Calculation time: seconds.
- Spin system parameters
- Basis set
- Algorithmic options
- Spinach housekeeping
- Sequence parameters
- Simulation
- Apodisation
- F2 Fourier transform
- Form States signal

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `liquid()`, `apodisation()`, `fftshift()`, `conj()`, `kfigure()`, `scale_figure()`, `plot_2d()`.
