# examples/nmr_spen/psyche_strychnine.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_spen/psyche_strychnine.m`
- Signature: `psyche_strychnine()`
- Total lines: 101

## Purpose

PSYCHE pure-shift NMR spectrum of strychnine. Calculation time: hours, faster on a GPU.

## Physical / mathematical content

- SPEN / ultrafast NMR examples. These files encode spatially dependent phase evolution and acquisition, linking pulse gradients, diffusion attenuation, and single-scan multidimensional encoding.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 10-11: Strychnine spin system; implemented by `[sys,inter]=strychnine({'1H'})`.
- Lines 13-14: Move shifts into spectral window; implemented by `inter.zeeman.scalar=inter.zeeman.scalar-5.0`.
- Lines 16-17: Magnetic induction; implemented by `sys.magnet=11.7`.
- Lines 19-20: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 25-26: Algorithmic options; implemented by `sys.tols.inter_cutoff=2.0`.
- Lines 30-31: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 34-35: Sequence parameters; implemented by `parameters.offset=0`.
- Lines 43-44: Saltire chirp parameters; implemented by `parameters.beta=20`.
- Lines 51-52: Coherent evolution timesteps; implemented by `parameters.timestep1=1/parameters.sweep(1)`.
- Lines 56-57: Sample parameters; implemented by `parameters.dims=0.015`.
- Lines 61-62: Diffusion and flow; implemented by `parameters.diff=0`.
- Lines 65-66: Initial and detection state phantoms; implemented by `parameters.rho0_ph={ones(parameters.npts,1)}`.
- Lines 71-72: Relaxation phantom; implemented by `parameters.rlx_ph={}`.
- Lines 75-76: Simulation; implemented by `fid=imaging(spin_system,@psyche,parameters)`.
- Lines 78-79: Reconstructing the pure shift FID; implemented by `np_chunk=parameters.sweep(2)/parameters.sweep(1)`.
- Lines 82-83: Apodisation; implemented by `fidps=apodisation(spin_system,fidps,{{'gauss',6}})`.
- Lines 85-87: Fourier transform; implemented by `spectrum_2d=fftshift(fft2(fid,parameters.zerofill(2), parameters.zerofill(1)))`.
- Lines 90-91: Plotting: full 2D version; implemented by `kfigure(); scale_figure([2.0 1.5]); subplot(1,2,1)`.

### Key state/data transformations

- Lines 11: computes `[sys,inter]` using `[sys,inter]=strychnine({'1H'})`.
- Lines 14: computes `inter.zeeman.scalar` using `inter.zeeman.scalar=inter.zeeman.scalar-5.0`.
- Lines 17: computes `sys.magnet` using `sys.magnet=11.7`.
- Lines 20: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 21: computes `bas.connectivity` using `bas.connectivity='scalar_couplings'`.
- Lines 22: computes `bas.approximation` using `bas.approximation='IK-2'`.
- Lines 23: computes `bas.space_level` using `bas.space_level=1`.
- Lines 26: computes `sys.tols.inter_cutoff` using `sys.tols.inter_cutoff=2.0`.
- Lines 27: computes `sys.disable` using `sys.disable={'pt','colorbar'}`.
- Lines 28: computes `sys.enable` using `sys.enable={'greedy'}`.
- Lines 31: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 35: computes `parameters.offset` using `parameters.offset=0`.
- Lines 36: computes `parameters.sweep` using `parameters.sweep=[100 5000]`.
- Lines 37: computes `parameters.npoints` using `parameters.npoints=[32 2048]`.
- Lines 38: computes `parameters.zerofill` using `parameters.zerofill=[128 8192]`.
- Lines 39: computes `parameters.spins` using `parameters.spins={'1H'}`.
- Lines 40: computes `parameters.axis_units` using `parameters.axis_units='Hz'`.
- Lines 41: computes `parameters.g_amp` using `parameters.g_amp=0.015`.

## Implementation structure

- PSYCHE pure-shift NMR spectrum of strychnine.
- Calculation time: hours, faster on a GPU.
- Strychnine spin system
- Move shifts into spectral window
- Magnetic induction
- Basis set
- Algorithmic options
- Spinach housekeeping
- Sequence parameters
- Saltire chirp parameters
- Coherent evolution timesteps
- Sample parameters

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `strychnine()`, `create()`, `basis()`, `state()`, `imaging()`, `fid()`, `fidps()`, `apodisation()`, `fftshift()`, `fft2()`, `kfigure()`, `scale_figure()`, `subplot()`, `plot_2d()`, `plot_1d()`.
