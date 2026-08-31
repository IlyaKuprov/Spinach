# examples/nmr_spen/psycosy_dbpa.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_spen/psycosy_dbpa.m`
- Signature: `psycosy_dbpa()`
- Total lines: 90

## Purpose

PSYCOSY of DBPA (dibromopropionic acid) ring. Calculation time: minutes on NVidia Tesla A100, much longer on CPU

## Physical / mathematical content

- SPEN / ultrafast NMR examples. These files encode spatially dependent phase evolution and acquisition, linking pulse gradients, diffusion attenuation, and single-scan multidimensional encoding.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 9-10: Magnet; implemented by `sys.magnet=14.1`.
- Lines 12-13: Spin system; implemented by `sys.isotopes={'1H','1H','1H', '1H'}`.
- Lines 15-16: Interactions; implemented by `inter.zeeman.scalar={4.49 3.9 3.7 4.2}`.
- Lines 23-24: Algorithmic options; implemented by `sys.enable={'greedy'}`.
- Lines 28-29: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 34-35: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 38-39: Sample geometry; implemented by `parameters.dims=0.015`.
- Lines 43-44: Diffusion and flow; implemented by `parameters.u=zeros(parameters.npts,1)`.
- Lines 47-48: Relaxation phantom; implemented by `parameters.rlx_ph={zeros(parameters.npts,1)}`.
- Lines 51-52: Initial and detection state phantoms; implemented by `parameters.rho0_ph={ones(parameters.npts,1)}`.
- Lines 57-58: Sequence parameters; implemented by `parameters.offset=2460`.
- Lines 67-68: Saltire chirp parameters; implemented by `parameters.sal_ang=20`.
- Lines 75-76: Simulation; implemented by `fid=imaging(spin_system,@psycosy,parameters)`.
- Lines 78-79: Apodisation; implemented by `fid=apodisation(spin_system,fid,{{'sqsin'},{'sqsin'}})`.
- Lines 81-82: Fourier transform; implemented by `spectrum=fftshift(fft2(fid,parameters.zerofill(2),parameters.zerofill(1)))`.
- Lines 84-85: Plotting and saving; implemented by `kfigure(); scale_figure([1.5 2.0])`.

### Key state/data transformations

- Lines 10: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 13: computes `sys.isotopes` using `sys.isotopes={'1H','1H','1H', '1H'}`.
- Lines 16: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={4.49 3.9 3.7 4.2}`.
- Lines 17: computes `inter.coupling.scalar{1,2}` using `inter.coupling.scalar{1,2}=11.3`.
- Lines 18: computes `inter.coupling.scalar{2,3}` using `inter.coupling.scalar{2,3}=10.1`.
- Lines 19: computes `inter.coupling.scalar{1,3}` using `inter.coupling.scalar{1,3}=4.3`.
- Lines 20: computes `inter.coupling.scalar{4,2}` using `inter.coupling.scalar{4,2}=0`.
- Lines 21: computes `inter.coupling.scalar{4,4}` using `inter.coupling.scalar{4,4}=0`.
- Lines 24: computes `sys.enable` using `sys.enable={'greedy'}`.
- Lines 25: computes `sys.tols.prox_cutoff` using `sys.tols.prox_cutoff=4.0`.
- Lines 26: computes `sys.tols.merge_dim` using `sys.tols.merge_dim=500`.
- Lines 29: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 30: computes `bas.approximation` using `bas.approximation='IK-2'`.
- Lines 31: computes `bas.space_level` using `bas.space_level=1`.
- Lines 32: computes `bas.connectivity` using `bas.connectivity='scalar_couplings'`.
- Lines 35: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 39: computes `parameters.dims` using `parameters.dims=0.015`.
- Lines 40: computes `parameters.npts` using `parameters.npts=100`.

## Implementation structure

- PSYCOSY of DBPA (dibromopropionic acid) ring.
- Calculation time: minutes on NVidia Tesla A100, much longer on CPU
- Magnet
- Spin system
- Interactions
- Algorithmic options
- Basis set
- Spinach housekeeping
- Sample geometry
- Diffusion and flow
- Relaxation phantom
- Initial and detection state phantoms

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `relaxation()`, `state()`, `imaging()`, `apodisation()`, `fftshift()`, `fft2()`, `kfigure()`, `scale_figure()`, `plot_2d()`.
