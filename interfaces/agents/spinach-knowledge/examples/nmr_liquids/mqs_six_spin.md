# examples/nmr_liquids/mqs_six_spin.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_liquids/mqs_six_spin.m`
- Signature: `mqs_six_spin()`
- Total lines: 90

## Purpose

Multiple-quantum (MQ) NMR experiment for a coupled system of six spins. Calculation time: minutes

## Physical / mathematical content

- Liquid-state NMR examples. The physics is scalar-coupling-mediated coherence transfer in weakly or moderately coupled spin systems, often in Liouville space. Typical mechanisms include INEPT-style polarisation transfer, J-refocusing, phase cycling, indirect evolution, and multidimensional detection.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 11-12: Magnetic field; implemented by `sys.magnet=14.1`.
- Lines 14-15: Chemical shifts; implemented by `sys.isotopes={'1H','1H','1H','1H','1H','1H'}`.
- Lines 18-19: 3J couplings; implemented by `inter.coupling.scalar{1,2}=8.77`.
- Lines 25-26: 4J couplings; implemented by `inter.coupling.scalar{1,3}=2.60`.
- Lines 31-32: 5J couplings; implemented by `inter.coupling.scalar{1,4}=2.00`.
- Lines 37-38: 6J couplings; implemented by `inter.coupling.scalar{1,5}=2.00`.
- Lines 42-43: Coherence to select; implemented by `parameters.mqorder=[+6 -1]`.
- Lines 45-46: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 49-50: Algorithmic options; implemented by `sys.enable={'greedy'}`.
- Lines 52-53: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 56-57: Initial and detection states; implemented by `parameters.rho0=state(spin_system,'Lz','1H')`.
- Lines 60-61: Sequence parameters; implemented by `parameters.angle=pi/2`.
- Lines 69-70: Tau delays; implemented by `parameters.delay_1=0.1`.
- Lines 73-74: Run simulation; implemented by `fid=liquid(spin_system,@mqs_refocus,parameters,'nmr')`.
- Lines 76-77: Apodisation; implemented by `fid=apodisation(spin_system,fid,{{'sqcos'},{'sqcos'}})`.
- Lines 79-81: Fourier transform; implemented by `spectrum_conv=fftshift(fft2(fid,parameters.zerofill(2), parameters.zerofill(1)))`.
- Lines 83-84: Plotting; implemented by `kfigure(); scale_figure([1.5 2.0])`.

### Key state/data transformations

- Lines 12: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 15: computes `sys.isotopes` using `sys.isotopes={'1H','1H','1H','1H','1H','1H'}`.
- Lines 16: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={-1.0 -0.5 0.0 +0.3 +0.7 +1.1}`.
- Lines 19: computes `inter.coupling.scalar{1,2}` using `inter.coupling.scalar{1,2}=8.77`.
- Lines 20: computes `inter.coupling.scalar{2,3}` using `inter.coupling.scalar{2,3}=7.40`.
- Lines 21: computes `inter.coupling.scalar{3,4}` using `inter.coupling.scalar{3,4}=7.40`.
- Lines 22: computes `inter.coupling.scalar{4,5}` using `inter.coupling.scalar{4,5}=8.77`.
- Lines 23: computes `inter.coupling.scalar{5,6}` using `inter.coupling.scalar{5,6}=8.77`.
- Lines 26: computes `inter.coupling.scalar{1,3}` using `inter.coupling.scalar{1,3}=2.60`.
- Lines 27: computes `inter.coupling.scalar{2,4}` using `inter.coupling.scalar{2,4}=2.60`.
- Lines 28: computes `inter.coupling.scalar{3,5}` using `inter.coupling.scalar{3,5}=2.60`.
- Lines 29: computes `inter.coupling.scalar{1,6}` using `inter.coupling.scalar{1,6}=2.60`.
- Lines 32: computes `inter.coupling.scalar{1,4}` using `inter.coupling.scalar{1,4}=2.00`.
- Lines 33: computes `inter.coupling.scalar{2,5}` using `inter.coupling.scalar{2,5}=2.00`.
- Lines 34: computes `inter.coupling.scalar{2,6}` using `inter.coupling.scalar{2,6}=2.00`.
- Lines 35: computes `inter.coupling.scalar{4,6}` using `inter.coupling.scalar{4,6}=2.00`.
- Lines 38: computes `inter.coupling.scalar{1,5}` using `inter.coupling.scalar{1,5}=2.00`.
- Lines 39: computes `inter.coupling.scalar{3,6}` using `inter.coupling.scalar{3,6}=2.00`.

## Implementation structure

- Multiple-quantum (MQ) NMR experiment for a coupled system
- of six spins.
- Calculation time: minutes
- Magnetic field
- Chemical shifts
- 3J couplings
- 4J couplings
- 5J couplings
- 6J couplings
- Coherence to select
- Basis set
- Algorithmic options

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `liquid()`, `apodisation()`, `fftshift()`, `fft2()`, `kfigure()`, `scale_figure()`, `plot_2d()`, `kxlabel()`, `kylabel()`.
