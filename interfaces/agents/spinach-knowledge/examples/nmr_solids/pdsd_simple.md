# examples/nmr_solids/pdsd_simple.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_solids/pdsd_simple.m`
- Signature: `pdsd_simple()`
- Total lines: 74

## Purpose

13C 2D PDSD spectrum of a simple test spin system. Calculation time: minutes, much faster on GPU.

## Physical / mathematical content

- Solid-state NMR examples. The key physics is anisotropic spin interactions under static or magic-angle-spinning conditions: chemical-shift anisotropy, dipolar coupling, quadrupolar coupling, cross-polarisation, and orientation averaging using Floquet, Fokker-Planck, or direct powder quadrature formalisms.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 10-11: Magnet field; implemented by `sys.magnet=21.1356`.
- Lines 13-14: Isotopes; implemented by `sys.isotopes={'13C','1H','1H','13C'}`.
- Lines 16-17: Interactions (HCCH fragment); implemented by `inter.zeeman.scalar={20.0 5.0 2.0 35.0}`.
- Lines 23-24: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 27-28: Algorithmic options; implemented by `sys.enable={'prop_cache'}`.
- Lines 30-31: Create the spin system structure; implemented by `spin_system=create(sys,inter)`.
- Lines 33-34: Build the basis; implemented by `spin_system=basis(spin_system,bas)`.
- Lines 36-37: Experiment parameters; implemented by `parameters.spins={'1H','13C'}`.
- Lines 48-49: Simulation; implemented by `fid=singlerot(spin_system,@pdsd,parameters,'nmr')`.
- Lines 51-52: Apodisation; implemented by `fid.cos=apodisation(spin_system,fid.cos,{{'sqcos'},{'sqcos'}})`.
- Lines 55-56: F2 Fourier transform; implemented by `f1_cos=imag(fftshift(fft(fid.cos,parameters.zerofill(2),1),1))`.
- Lines 59-60: States signal; implemented by `f1_states=f1_cos-1i*f1_sin`.
- Lines 62-63: F1 Fourier transform and real part; implemented by `spectrum=fftshift(fft(f1_states,parameters.zerofill(1),2),2)`.
- Lines 66-67: Plotting; implemented by `kfigure(); scale_figure([1.5 2.0])`.

### Key state/data transformations

- Lines 11: computes `sys.magnet` using `sys.magnet=21.1356`.
- Lines 14: computes `sys.isotopes` using `sys.isotopes={'13C','1H','1H','13C'}`.
- Lines 17: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={20.0 5.0 2.0 35.0}`.
- Lines 18-21: computes `inter.coordinates` using `inter.coordinates={[-2.26 0.15 0.00], [-1.90 0.66 -0.87], [-0.67 0.88 1.26], [-1.74 0.88 1.26]}`.
- Lines 24: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 25: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 28: computes `sys.enable` using `sys.enable={'prop_cache'}`.
- Lines 31: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 37: computes `parameters.spins` using `parameters.spins={'1H','13C'}`.
- Lines 38: computes `parameters.rate` using `parameters.rate=10000`.
- Lines 39: computes `parameters.tmix` using `parameters.tmix=10e-3`.
- Lines 40: computes `parameters.max_rank` using `parameters.max_rank=11`.
- Lines 41: computes `parameters.axis` using `parameters.axis=[sqrt(2/3) 0 sqrt(1/3)]`.
- Lines 42: computes `parameters.grid` using `parameters.grid='rep_2ang_100pts_sph'`.
- Lines 43: computes `parameters.npoints` using `parameters.npoints=[256 256]`.
- Lines 44: computes `parameters.zerofill` using `parameters.zerofill=[1024 1024]`.
- Lines 45: computes `parameters.offset` using `parameters.offset=[3150 6750]`.
- Lines 46: computes `parameters.sweep` using `parameters.sweep=10000`.

## Implementation structure

- 13C 2D PDSD spectrum of a simple test spin system.
- Calculation time: minutes, much faster on GPU.
- Magnet field
- Isotopes
- Interactions (HCCH fragment)
- Basis set
- Algorithmic options
- Create the spin system structure
- Build the basis
- Experiment parameters
- Simulation
- Apodisation

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `singlerot()`, `apodisation()`, `fftshift()`, `kfigure()`, `scale_figure()`, `plot_2d()`.
