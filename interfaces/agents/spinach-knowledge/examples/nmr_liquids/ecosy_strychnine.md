# examples/nmr_liquids/ecosy_strychnine.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_liquids/ecosy_strychnine.m`
- Signature: `ecosy_strychnine()`
- Total lines: 60

## Purpose

E.COSY spectrum of strychnine. Calculation time: minutes

## Physical / mathematical content

- Liquid-state NMR examples. The physics is scalar-coupling-mediated coherence transfer in weakly or moderately coupled spin systems, often in Liouville space. Typical mechanisms include INEPT-style polarisation transfer, J-refocusing, phase cycling, indirect evolution, and multidimensional detection.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 9-10: Read the spin system properties; implemented by `[sys,inter]=strychnine({'1H'})`.
- Lines 12-13: Magnet field; implemented by `sys.magnet=5.9`.
- Lines 15-16: Algorithmic options; implemented by `sys.enable={'greedy'}`.
- Lines 19-20: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 25-26: Sequence parameters; implemented by `parameters.offset=1200`.
- Lines 33-34: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 37-38: Simulation; implemented by `fid=liquid(spin_system,@ecosy,parameters,'nmr')`.
- Lines 40-41: Apodisation; implemented by `fid.cos=apodisation(spin_system,fid.cos,{{'sqcos'},{'sqcos'}})`.
- Lines 44-45: F2 Fourier transform; implemented by `f1_cos=fftshift(fft(fid.cos,parameters.zerofill(2),1),1)`.
- Lines 48-49: Form States signal; implemented by `f1_states=real(f1_cos)+1i*imag(f1_sin)`.
- Lines 51-52: F1 Fourier transform; implemented by `spectrum=fftshift(fft(f1_states,parameters.zerofill(1),2),2)`.
- Lines 54-55: Plotting; implemented by `kfigure(); scale_figure([1.5 1.5])`.

### Key state/data transformations

- Lines 10: computes `[sys,inter]` using `[sys,inter]=strychnine({'1H'})`.
- Lines 13: computes `sys.magnet` using `sys.magnet=5.9`.
- Lines 16: computes `sys.enable` using `sys.enable={'greedy'}`.
- Lines 17: computes `sys.tols.prox_cutoff` using `sys.tols.prox_cutoff=4.0`.
- Lines 20: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 21: computes `bas.approximation` using `bas.approximation='IK-2'`.
- Lines 22: computes `bas.connectivity` using `bas.connectivity='scalar_couplings'`.
- Lines 23: computes `bas.space_level` using `bas.space_level=1`.
- Lines 26: computes `parameters.offset` using `parameters.offset=1200`.
- Lines 27: computes `parameters.sweep` using `parameters.sweep=2200`.
- Lines 28: computes `parameters.npoints` using `parameters.npoints=[512 512]`.
- Lines 29: computes `parameters.zerofill` using `parameters.zerofill=[2048 2048]`.
- Lines 30: computes `parameters.spins` using `parameters.spins={'1H'}`.
- Lines 31: computes `parameters.axis_units` using `parameters.axis_units='ppm'`.
- Lines 34: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 38: computes `fid` using `fid=liquid(spin_system,@ecosy,parameters,'nmr')`.
- Lines 41: computes `fid.cos` using `fid.cos=apodisation(spin_system,fid.cos,{{'sqcos'},{'sqcos'}})`.
- Lines 42: computes `fid.sin` using `fid.sin=apodisation(spin_system,fid.sin,{{'sqcos'},{'sqcos'}})`.

## Implementation structure

- E.COSY spectrum of strychnine.
- Calculation time: minutes
- Read the spin system properties
- Magnet field
- Algorithmic options
- Basis set
- Sequence parameters
- Spinach housekeeping
- Simulation
- Apodisation
- F2 Fourier transform
- Form States signal

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `strychnine()`, `create()`, `basis()`, `liquid()`, `apodisation()`, `fftshift()`, `kfigure()`, `scale_figure()`, `plot_2d()`.
