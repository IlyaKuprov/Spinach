# examples/nmr_spen/ufcosy_3spin.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_spen/ufcosy_3spin.m`
- Signature: `ufcosy_3spin()`
- Total lines: 77

## Purpose

Ultrafast COSY for a coupled three-spin system. Calculation time: minutes on NVidia Tesla A100, much longer on CPU Jean-Nicolas Dumez Ludmilla Guduff

## Physical / mathematical content

- SPEN / ultrafast NMR examples. These files encode spatially dependent phase evolution and acquisition, linking pulse gradients, diffusion attenuation, and single-scan multidimensional encoding.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 10-11: Interactions; implemented by `sys.magnet=14.0`.
- Lines 19-20: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 23-24: Algorithmic options; implemented by `sys.disable={'pt'}`.
- Lines 27-28: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 31-32: Sample geometry; implemented by `parameters.dims=0.015`.
- Lines 36-37: Relaxation phantom; implemented by `parameters.rlx_ph={}`.
- Lines 40-41: Initial and detection state phantoms; implemented by `parameters.rho0_ph={ones(parameters.npts,1)}`.
- Lines 46-47: Diffusion and flow; implemented by `parameters.u=zeros(parameters.npts,1)`.
- Lines 50-51: Acquisition parameters; implemented by `parameters.spins={'1H'}`.
- Lines 58-59: Encoding parameters; implemented by `parameters.pulsenpoints=1000`.
- Lines 65-66: Coherence selection parameters; implemented by `parameters.Gp=0.47`.
- Lines 69-70: Simulation; implemented by `fid=imaging(spin_system,@spencosy,parameters)`.
- Lines 72-73: Processing and plotting; implemented by `spectrum=fftshift(fft(fid,[],2),2)`.

### Key state/data transformations

- Lines 11: computes `sys.magnet` using `sys.magnet=14.0`.
- Lines 12: computes `sys.isotopes` using `sys.isotopes={'1H','1H','1H'}`.
- Lines 13: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={3.70 3.92 4.50}-{4.0 4.0 4.0}`.
- Lines 14: computes `inter.coupling.scalar{1,2}` using `inter.coupling.scalar{1,2}=10`.
- Lines 15: computes `inter.coupling.scalar{2,3}` using `inter.coupling.scalar{2,3}=12`.
- Lines 16: computes `inter.coupling.scalar{1,3}` using `inter.coupling.scalar{1,3}=4`.
- Lines 17: computes `inter.coupling.scalar{3,3}` using `inter.coupling.scalar{3,3}=0`.
- Lines 20: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 21: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 24: computes `sys.disable` using `sys.disable={'pt'}`.
- Lines 25: computes `sys.enable` using `sys.enable={'greedy'}`.
- Lines 28: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 32: computes `parameters.dims` using `parameters.dims=0.015`.
- Lines 33: computes `parameters.npts` using `parameters.npts=500`.
- Lines 34: computes `parameters.deriv` using `parameters.deriv={'period',3}`.
- Lines 37: computes `parameters.rlx_ph` using `parameters.rlx_ph={}`.
- Lines 38: computes `parameters.rlx_op` using `parameters.rlx_op={}`.
- Lines 41: computes `parameters.rho0_ph` using `parameters.rho0_ph={ones(parameters.npts,1)}`.

## Implementation structure

- Ultrafast COSY for a coupled three-spin system.
- Calculation time: minutes on NVidia Tesla A100, much longer on CPU
- Jean-Nicolas Dumez
- Ludmilla Guduff
- Interactions
- Basis set
- Algorithmic options
- Spinach housekeeping
- Sample geometry
- Relaxation phantom
- Initial and detection state phantoms
- Diffusion and flow

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `imaging()`, `fftshift()`, `kfigure()`.
