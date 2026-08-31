# examples/nmr_liquids/gcosy90_strychnine.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_liquids/gcosy90_strychnine.m`
- Signature: `gcosy90_strychnine()`
- Total lines: 67

## Purpose

Gradient-selected COSY spectrum of strychnine. Calculation time: minutes

## Physical / mathematical content

- Liquid-state NMR examples. The physics is scalar-coupling-mediated coherence transfer in weakly or moderately coupled spin systems, often in Liouville space. Typical mechanisms include INEPT-style polarisation transfer, J-refocusing, phase cycling, indirect evolution, and multidimensional detection.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 11-12: Read the spin system properties; implemented by `[sys,inter]=strychnine({'1H'})`.
- Lines 14-15: Magnet field; implemented by `sys.magnet=5.9`.
- Lines 17-18: Algorithmic options; implemented by `sys.enable={'greedy'}`.
- Lines 21-22: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 27-28: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 31-32: Sequence parameters; implemented by `parameters.angle=pi/2`.
- Lines 45-46: Simulation; implemented by `fid=liquid(spin_system,@gcosy,parameters,'nmr')`.
- Lines 48-49: Apodisation; implemented by `fid.pos=apodisation(spin_system,fid.pos,{{'sqcos'},{'sqcos'}})`.
- Lines 52-53: F2 Fourier transform; implemented by `f1_pos=fftshift(fft(fid.pos,parameters.zerofill(2),1),1)`.
- Lines 56-57: Form echo/anti-echo signal; implemented by `fid=f1_pos+conj(f1_neg)`.
- Lines 59-60: F1 Fourier transform; implemented by `spectrum=fftshift(fft(fid,parameters.zerofill(1),2),2)`.
- Lines 62-63: Plotting; implemented by `kfigure(); scale_figure([1.5 2.0])`.

### Key state/data transformations

- Lines 12: computes `[sys,inter]` using `[sys,inter]=strychnine({'1H'})`.
- Lines 15: computes `sys.magnet` using `sys.magnet=5.9`.
- Lines 18: computes `sys.enable` using `sys.enable={'greedy'}`.
- Lines 19: computes `sys.tols.prox_cutoff` using `sys.tols.prox_cutoff=4.0`.
- Lines 22: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 23: computes `bas.approximation` using `bas.approximation='IK-2'`.
- Lines 24: computes `bas.connectivity` using `bas.connectivity='scalar_couplings'`.
- Lines 25: computes `bas.space_level` using `bas.space_level=1`.
- Lines 28: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 32: computes `parameters.angle` using `parameters.angle=pi/2`.
- Lines 33: computes `parameters.offset` using `parameters.offset=1200`.
- Lines 34: computes `parameters.sweep` using `parameters.sweep=2200`.
- Lines 35: computes `parameters.npoints` using `parameters.npoints=[512 512]`.
- Lines 36: computes `parameters.zerofill` using `parameters.zerofill=[2048 2048]`.
- Lines 37: computes `parameters.spins` using `parameters.spins={'1H'}`.
- Lines 38: computes `parameters.axis_units` using `parameters.axis_units='ppm'`.
- Lines 39: computes `parameters.g_amp` using `parameters.g_amp=3`.
- Lines 40: computes `parameters.g_dur` using `parameters.g_dur=2e-3`.

## Implementation structure

- Gradient-selected COSY spectrum of strychnine.
- Calculation time: minutes
- Read the spin system properties
- Magnet field
- Algorithmic options
- Basis set
- Spinach housekeeping
- Sequence parameters
- Simulation
- Apodisation
- F2 Fourier transform
- Form echo/anti-echo signal

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `strychnine()`, `create()`, `basis()`, `liquid()`, `apodisation()`, `fftshift()`, `conj()`, `kfigure()`, `scale_figure()`, `plot_2d()`.
