# examples/nmr_liquids/mqs_propanol.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_liquids/mqs_propanol.m`
- Signature: `mqs_propanol()`
- Total lines: 101

## Purpose

Multiple-quantum NMR experiment for a propanol spin system. Calculation time: seconds

## Physical / mathematical content

- Liquid-state NMR examples. The physics is scalar-coupling-mediated coherence transfer in weakly or moderately coupled spin systems, often in Liouville space. Typical mechanisms include INEPT-style polarisation transfer, J-refocusing, phase cycling, indirect evolution, and multidimensional detection.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 12-13: Magnet field; implemented by `sys.magnet=14.1`.
- Lines 15-16: Chemical shifts; implemented by `sys.isotopes={'1H','1H','1H','1H','1H','1H','1H'}`.
- Lines 19-20: 2J couplings; implemented by `inter.coupling.scalar{1,3}=7.5`.
- Lines 31-32: 3J couplings; implemented by `inter.coupling.scalar{1,5}=0.5`.
- Lines 40-41: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 46-47: Algorithmic options; implemented by `sys.enable={'greedy'}`.
- Lines 49-50: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 53-54: Initial and detection states; implemented by `parameters.rho0=state(spin_system,'Lz','1H')`.
- Lines 57-58: Sequence parameters; implemented by `parameters.angle=pi/2`.
- Lines 66-67: Tau delays; implemented by `tau_max=[0.0333 0.0710 0.5000]`.
- Lines 69-70: Coherence to select; implemented by `parameters.mqorder=+3`.
- Lines 72-73: Get a figure going; implemented by `kfigure(); scale_figure([2.5 1.5])`.
- Lines 75-76: Loop over tau delays; implemented by `for n=1:numel(tau_max)`.
- Lines 78-79: tau delay; implemented by `parameters.delay=tau_max(n)`.
- Lines 81-82: Run simulation; implemented by `fid=liquid(spin_system,@mqs,parameters,'nmr')`.
- Lines 84-85: Apodisation; implemented by `fid=apodisation(spin_system,fid,{{'sin'},{'sin'}})`.
- Lines 87-89: Fourier transform; implemented by `spectrum=fftshift(fft2(fid,parameters.zerofill(2), parameters.zerofill(1)))`.
- Lines 91-92: Get subplots for differnt tau values; implemented by `subplot(1,3,n)`.

### Control flow inferred from the code

- Line 76: `for` loop over `n=1:numel(tau_max)`.

### Key state/data transformations

- Lines 13: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 16: computes `sys.isotopes` using `sys.isotopes={'1H','1H','1H','1H','1H','1H','1H'}`.
- Lines 17: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={3.438,3.438,1.429,1.429,0.775,0.775,0.775}`.
- Lines 20: computes `inter.coupling.scalar{1,3}` using `inter.coupling.scalar{1,3}=7.5`.
- Lines 21: computes `inter.coupling.scalar{1,4}` using `inter.coupling.scalar{1,4}=7.5`.
- Lines 22: computes `inter.coupling.scalar{2,3}` using `inter.coupling.scalar{2,3}=7.5`.
- Lines 23: computes `inter.coupling.scalar{2,4}` using `inter.coupling.scalar{2,4}=7.5`.
- Lines 24: computes `inter.coupling.scalar{3,5}` using `inter.coupling.scalar{3,5}=7.0`.
- Lines 25: computes `inter.coupling.scalar{3,6}` using `inter.coupling.scalar{3,6}=7.0`.
- Lines 26: computes `inter.coupling.scalar{3,7}` using `inter.coupling.scalar{3,7}=7.0`.
- Lines 27: computes `inter.coupling.scalar{4,5}` using `inter.coupling.scalar{4,5}=7.0`.
- Lines 28: computes `inter.coupling.scalar{4,6}` using `inter.coupling.scalar{4,6}=7.0`.
- Lines 29: computes `inter.coupling.scalar{4,7}` using `inter.coupling.scalar{4,7}=7.0`.
- Lines 32: computes `inter.coupling.scalar{1,5}` using `inter.coupling.scalar{1,5}=0.5`.
- Lines 33: computes `inter.coupling.scalar{1,6}` using `inter.coupling.scalar{1,6}=0.5`.
- Lines 34: computes `inter.coupling.scalar{1,7}` using `inter.coupling.scalar{1,7}=0.5`.
- Lines 35: computes `inter.coupling.scalar{2,5}` using `inter.coupling.scalar{2,5}=0.5`.
- Lines 36: computes `inter.coupling.scalar{2,6}` using `inter.coupling.scalar{2,6}=0.5`.

## Implementation structure

- Multiple-quantum NMR experiment for a propanol
- spin system.
- Calculation time: seconds
- Magnet field
- Chemical shifts
- 2J couplings
- 3J couplings
- Basis set
- Algorithmic options
- Spinach housekeeping
- Initial and detection states
- Sequence parameters

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `kfigure()`, `scale_figure()`, `tau_max()`, `liquid()`, `apodisation()`, `fftshift()`, `fft2()`, `subplot()`, `plot_2d()`, `num2str()`.
