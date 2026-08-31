# examples/nmr_liquids/dqf_cosy_sucrose.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_liquids/dqf_cosy_sucrose.m`
- Signature: `dqf_cosy_sucrose()`
- Total lines: 62

## Purpose

DQF-COSY spectrum of sucrose (magnetic parameters computed with DFT). Calculation time: minutes

## Physical / mathematical content

- Liquid-state NMR examples. The physics is scalar-coupling-mediated coherence transfer in weakly or moderately coupled spin systems, often in Liouville space. Typical mechanisms include INEPT-style polarisation transfer, J-refocusing, phase cycling, indirect evolution, and multidimensional detection.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 10-11: Spin system properties (vacuum DFT calculation); implemented by `options.min_j=2.0; options.no_xyz=1`.
- Lines 14-15: Magnet field; implemented by `sys.magnet=5.9`.
- Lines 17-18: Algorithmic options; implemented by `sys.enable={'greedy'}`.
- Lines 21-22: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 27-28: Sequence parameters; implemented by `parameters.offset=800`.
- Lines 35-36: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 39-40: Simulation; implemented by `fid=liquid(spin_system,@dqf_cosy,parameters,'nmr')`.
- Lines 42-43: Apodization; implemented by `fid.cos=apodisation(spin_system,fid.cos,{{'cos'},{'cos'}})`.
- Lines 46-47: F2 Fourier transform; implemented by `f1_cos=fftshift(fft(fid.cos,parameters.zerofill(2),1),1)`.
- Lines 50-51: Form States signal; implemented by `f1_states=real(f1_cos)-1i*real(f1_sin)`.
- Lines 53-54: F1 Fourier transform; implemented by `spectrum=fftshift(fft(f1_states,parameters.zerofill(1),2),2)`.
- Lines 56-57: Plotting; implemented by `kfigure(); scale_figure([1.5 2.0])`.

### Key state/data transformations

- Lines 11: computes `options.min_j` using `options.min_j=2.0; options.no_xyz=1`.
- Lines 12-13: computes `[sys,inter]` using `[sys,inter]=g2spinach(gparse('../standard_systems/sucrose.log'), {{'H','1H'}},31.8,options)`.
- Lines 15: computes `sys.magnet` using `sys.magnet=5.9`.
- Lines 18: computes `sys.enable` using `sys.enable={'greedy'}`.
- Lines 19: computes `sys.tols.prox_cutoff` using `sys.tols.prox_cutoff=4.0`.
- Lines 22: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 23: computes `bas.approximation` using `bas.approximation='IK-2'`.
- Lines 24: computes `bas.connectivity` using `bas.connectivity='scalar_couplings'`.
- Lines 25: computes `bas.space_level` using `bas.space_level=1`.
- Lines 28: computes `parameters.offset` using `parameters.offset=800`.
- Lines 29: computes `parameters.sweep` using `parameters.sweep=1700`.
- Lines 30: computes `parameters.npoints` using `parameters.npoints=[512 512]`.
- Lines 31: computes `parameters.zerofill` using `parameters.zerofill=[2048 2048]`.
- Lines 32: computes `parameters.spins` using `parameters.spins={'1H'}`.
- Lines 33: computes `parameters.axis_units` using `parameters.axis_units='ppm'`.
- Lines 36: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 40: computes `fid` using `fid=liquid(spin_system,@dqf_cosy,parameters,'nmr')`.
- Lines 43: computes `fid.cos` using `fid.cos=apodisation(spin_system,fid.cos,{{'cos'},{'cos'}})`.

## Implementation structure

- DQF-COSY spectrum of sucrose (magnetic parameters computed with DFT).
- Calculation time: minutes
- Spin system properties (vacuum DFT calculation)
- Magnet field
- Algorithmic options
- Basis set
- Sequence parameters
- Spinach housekeeping
- Simulation
- Apodization
- F2 Fourier transform
- Form States signal

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `g2spinach()`, `gparse()`, `create()`, `basis()`, `liquid()`, `apodisation()`, `fftshift()`, `kfigure()`, `scale_figure()`, `plot_2d()`.
