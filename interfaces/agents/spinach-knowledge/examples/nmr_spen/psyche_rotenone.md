# examples/nmr_spen/psyche_rotenone.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_spen/psyche_rotenone.m`
- Signature: `psyche_rotenone()`
- Total lines: 123

## Purpose

PSYCHE pure-shift NMR spectrum of rotenone. Calculation time: hours, faster on a GPU.

## Physical / mathematical content

- SPEN / ultrafast NMR examples. These files encode spatially dependent phase evolution and acquisition, linking pulse gradients, diffusion attenuation, and single-scan multidimensional encoding.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 10-11: Magnetic induction; implemented by `sys.magnet=11.7`.
- Lines 13-16: Spin system; implemented by `sys.isotopes={'1H','1H','1H','1H','1H','1H','1H','1H','1H', '1H','1H','1H','1H','1H','1H','1H','1H','1H', '1H','1H','1H','1H'}`.
- Lines 18-21: Interactions; implemented by `inter.zeeman.scalar={6.72 6.40 4.13 4.56 4.89 6.46 7.79 3.79 2.91 3.27 5.19 4.89 5.03 1.72 1.72 1.72 3.72 3.72 3.72 3.76 3.76 3.76}-5.0`.
- Lines 42-43: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 48-49: Algorithmic options; implemented by `sys.disable={'pt','colorbar'}`.
- Lines 52-53: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 56-57: Sequence parameters; implemented by `parameters.offset=0`.
- Lines 65-66: Saltire chirp parameters; implemented by `parameters.beta=20`.
- Lines 73-74: Coherent evolution timesteps; implemented by `parameters.timestep1=1/parameters.sweep(1)`.
- Lines 78-79: Sample parameters; implemented by `parameters.dims=0.015`.
- Lines 83-84: Diffusion and flow; implemented by `parameters.diff=0`.
- Lines 87-88: Initial and detection state phantoms; implemented by `parameters.rho0_ph={ones(parameters.npts,1)}`.
- Lines 93-94: Relaxation phantom; implemented by `parameters.rlx_ph={}`.
- Lines 97-98: Simulation; implemented by `fid=imaging(spin_system,@psyche,parameters)`.
- Lines 100-101: Reconstructing the pure shift FID; implemented by `np_chunk=parameters.sweep(2)/parameters.sweep(1)`.
- Lines 104-105: Apodization; implemented by `fidps=apodisation(spin_system,fidps,{{'gauss',6}})`.
- Lines 107-109: Fourier transform; implemented by `spectrum_2d=fftshift(fft2(fid,parameters.zerofill(2), parameters.zerofill(1)))`.
- Lines 112-113: Plotting: full 2D version; implemented by `kfigure(); scale_figure([2.0 1.5]); subplot(1,2,1)`.

### Key state/data transformations

- Lines 11: computes `sys.magnet` using `sys.magnet=11.7`.
- Lines 14-16: computes `sys.isotopes` using `sys.isotopes={'1H','1H','1H','1H','1H','1H','1H','1H','1H', '1H','1H','1H','1H','1H','1H','1H','1H','1H', '1H','1H','1H','1H'}`.
- Lines 19-21: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={6.72 6.40 4.13 4.56 4.89 6.46 7.79 3.79 2.91 3.27 5.19 4.89 5.03 1.72 1.72 1.72 3.72 3.72 3.72 3.76 3.76 3.76}-5.0`.
- Lines 22: computes `inter.coupling.scalar{3,4}` using `inter.coupling.scalar{3,4}=12.1`.
- Lines 23: computes `inter.coupling.scalar{4,5}` using `inter.coupling.scalar{4,5}=3.1`.
- Lines 24: computes `inter.coupling.scalar{3,5}` using `inter.coupling.scalar{3,5}=1.0`.
- Lines 25: computes `inter.coupling.scalar{3,8}` using `inter.coupling.scalar{3,8}=1.0`.
- Lines 26: computes `inter.coupling.scalar{1,8}` using `inter.coupling.scalar{1,8}=1.0`.
- Lines 27: computes `inter.coupling.scalar{6,7}` using `inter.coupling.scalar{6,7}=8.6`.
- Lines 28: computes `inter.coupling.scalar{5,8}` using `inter.coupling.scalar{5,8}=4.1`.
- Lines 29: computes `inter.coupling.scalar{7,9}` using `inter.coupling.scalar{7,9}=0.7`.
- Lines 30: computes `inter.coupling.scalar{7,10}` using `inter.coupling.scalar{7,10}=0.7`.
- Lines 31: computes `inter.coupling.scalar{9,10}` using `inter.coupling.scalar{9,10}=15.8`.
- Lines 32: computes `inter.coupling.scalar{10,11}` using `inter.coupling.scalar{10,11}=9.8`.
- Lines 33: computes `inter.coupling.scalar{9,11}` using `inter.coupling.scalar{9,11}=8.1`.
- Lines 34: computes `inter.coupling.scalar{13,14}` using `inter.coupling.scalar{13,14}=1.5`.
- Lines 35: computes `inter.coupling.scalar{12,14}` using `inter.coupling.scalar{12,14}=0.9`.
- Lines 36: computes `inter.coupling.scalar{13,15}` using `inter.coupling.scalar{13,15}=1.5`.

## Implementation structure

- PSYCHE pure-shift NMR spectrum of rotenone.
- Calculation time: hours, faster on a GPU.
- Magnetic induction
- Spin system
- Interactions
- Basis set
- Algorithmic options
- Spinach housekeeping
- Sequence parameters
- Saltire chirp parameters
- Coherent evolution timesteps
- Sample parameters

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `imaging()`, `fid()`, `fidps()`, `apodisation()`, `fftshift()`, `fft2()`, `kfigure()`, `scale_figure()`, `subplot()`, `plot_2d()`, `plot_1d()`.
