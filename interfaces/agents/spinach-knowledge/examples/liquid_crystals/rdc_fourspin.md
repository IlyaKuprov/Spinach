# examples/liquid_crystals/rdc_fourspin.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/liquid_crystals/rdc_fourspin.m`
- Signature: `rdc_fourspin()`
- Total lines: 70

## Purpose

CLIP-HSQC spectrum of a four-spin system in a liquid crystal with a user-specified order matrix. Calculation times: seconds.

## Physical / mathematical content

- Liquid-crystal examples. These scripts exploit partial ordering and Saupe-tensor physics, so anisotropic couplings survive orientational averaging and generate residual dipolar couplings or anisotropic transfer behaviour.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 10-11: Magnet field; implemented by `sys.magnet=5.9`.
- Lines 13-14: Spin system and interactions; implemented by `sys.isotopes={'1H','1H','13C','13C'}`.
- Lines 26-27: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 30-31: Algorithmic options; implemented by `sys.disable={'colorbar'}`.
- Lines 33-34: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 37-38: Sequence parameters; implemented by `parameters.sweep=[5000 5000]`.
- Lines 47-48: Simulation; implemented by `fid=liquid(spin_system,@clip_hsqc,parameters,'nmr')`.
- Lines 50-51: Apodisation; implemented by `fid.pos=apodisation(spin_system,fid.pos,{{'sqcos'},{'sqcos'}})`.
- Lines 54-55: F2 Fourier transform; implemented by `f1_pos=fftshift(fft(fid.pos,parameters.zerofill(2),1),1)`.
- Lines 58-59: Form States signal; implemented by `fid=f1_pos+conj(f1_neg)`.
- Lines 61-62: F1 Fourier transform; implemented by `spectrum=fftshift(fft(fid,parameters.zerofill(1),2),2)`.
- Lines 64-65: Plotting; implemented by `kfigure(); scale_figure([1.5 2.0])`.

### Key state/data transformations

- Lines 11: computes `sys.magnet` using `sys.magnet=5.9`.
- Lines 14: computes `sys.isotopes` using `sys.isotopes={'1H','1H','13C','13C'}`.
- Lines 15: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={1,10,80,50}`.
- Lines 16: computes `inter.coupling.scalar` using `inter.coupling.scalar={0 2.5 250 50`.
- Lines 20: computes `inter.coordinates` using `inter.coordinates={[1.0 0.0 0.0]`.
- Lines 24: computes `inter.order_matrix` using `inter.order_matrix={diag([1e-3 2e-3 -3e-3])}`.
- Lines 27: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 28: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 31: computes `sys.disable` using `sys.disable={'colorbar'}`.
- Lines 34: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 38: computes `parameters.sweep` using `parameters.sweep=[5000 5000]`.
- Lines 39: computes `parameters.offset` using `parameters.offset=[4250 1200]`.
- Lines 40: computes `parameters.npoints` using `parameters.npoints=[128 128]`.
- Lines 41: computes `parameters.zerofill` using `parameters.zerofill=[512 512]`.
- Lines 42: computes `parameters.spins` using `parameters.spins={'13C','1H'}`.
- Lines 43: computes `parameters.J` using `parameters.J=140`.
- Lines 44: computes `parameters.needs` using `parameters.needs={'rdc'}`.
- Lines 45: computes `parameters.axis_units` using `parameters.axis_units='ppm'`.

## Implementation structure

- CLIP-HSQC spectrum of a four-spin system in a liquid crystal
- with a user-specified order matrix.
- Calculation times: seconds.
- Magnet field
- Spin system and interactions
- Basis set
- Algorithmic options
- Spinach housekeeping
- Sequence parameters
- Simulation
- Apodisation
- F2 Fourier transform

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `liquid()`, `apodisation()`, `fftshift()`, `conj()`, `kfigure()`, `scale_figure()`, `plot_2d()`.
