# examples/nmr_liquids/tocsy_sucrose.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_liquids/tocsy_sucrose.m`
- Signature: `tocsy_sucrose()`
- Total lines: 65

## Purpose

TOCSY spectrum of sucrose (magnetic parameters computed with DFT). Calculation time: seconds

## Physical / mathematical content

- Liquid-state NMR examples. The physics is scalar-coupling-mediated coherence transfer in weakly or moderately coupled spin systems, often in Liouville space. Typical mechanisms include INEPT-style polarisation transfer, J-refocusing, phase cycling, indirect evolution, and multidimensional detection.
- Propagation is accelerated with a Krylov-subspace method, replacing direct matrix exponentiation by projection into a much smaller Arnoldi/Lanczos-type subspace.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.
- A Krylov-subspace or Arnoldi construction is used to avoid forming or exponentiating very large dense propagators directly.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 9-10: Spin system properties (vacuum DFT calculation); implemented by `options.min_j=1.0`.
- Lines 13-14: Magnet field; implemented by `sys.magnet=5.9`.
- Lines 16-17: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 22-23: Algorithmic options; implemented by `sys.enable={'greedy'}`.
- Lines 27-28: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 31-32: Sequence parameters; implemented by `parameters.tmix=0.100`.
- Lines 42-43: Simulation; implemented by `fid=liquid(spin_system,@tocsy,parameters,'nmr')`.
- Lines 45-46: Apodisation; implemented by `fid.cos=apodisation(spin_system,fid.cos,{{'sqcos'},{'sqcos'}})`.
- Lines 49-50: F2 Fourier transform; implemented by `f1_cos=imag(fftshift(fft(fid.cos,parameters.zerofill(2),1),1))`.
- Lines 53-54: States signal; implemented by `f1_states=f1_cos-1i*f1_sin`.
- Lines 56-57: F1 Fourier transform; implemented by `spectrum=fftshift(fft(f1_states,parameters.zerofill(1),2),2)`.
- Lines 59-60: Plotting; implemented by `kfigure(); scale_figure([1.5 2.0])`.

### Key state/data transformations

- Lines 10: computes `options.min_j` using `options.min_j=1.0`.
- Lines 11-12: computes `[sys,inter]` using `[sys,inter]=g2spinach(gparse('../standard_systems/sucrose.log'), {{'H','1H'}},31.8,options)`.
- Lines 14: computes `sys.magnet` using `sys.magnet=5.9`.
- Lines 17: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 18: computes `bas.approximation` using `bas.approximation='IK-2'`.
- Lines 19: computes `bas.connectivity` using `bas.connectivity='scalar_couplings'`.
- Lines 20: computes `bas.space_level` using `bas.space_level=1`.
- Lines 23: computes `sys.enable` using `sys.enable={'greedy'}`.
- Lines 24: computes `sys.disable` using `sys.disable={'krylov'}`.
- Lines 25: computes `sys.tols.prox_cutoff` using `sys.tols.prox_cutoff=4.0`.
- Lines 28: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 32: computes `parameters.tmix` using `parameters.tmix=0.100`.
- Lines 33: computes `parameters.lamp` using `parameters.lamp=1e4`.
- Lines 34: computes `parameters.offset` using `parameters.offset=800`.
- Lines 35: computes `parameters.sweep` using `parameters.sweep=[1700 1700]`.
- Lines 36: computes `parameters.npoints` using `parameters.npoints=[512 512]`.
- Lines 37: computes `parameters.zerofill` using `parameters.zerofill=[2048 2048]`.
- Lines 38: computes `parameters.spins` using `parameters.spins={'1H'}`.

## Implementation structure

- TOCSY spectrum of sucrose (magnetic parameters computed with DFT).
- Calculation time: seconds
- Spin system properties (vacuum DFT calculation)
- Magnet field
- Basis set
- Algorithmic options
- Spinach housekeeping
- Sequence parameters
- Simulation
- Apodisation
- F2 Fourier transform
- States signal

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `g2spinach()`, `gparse()`, `create()`, `basis()`, `state()`, `liquid()`, `apodisation()`, `fftshift()`, `kfigure()`, `scale_figure()`, `plot_2d()`.
