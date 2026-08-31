# examples/nmr_spen/ufdosycosy_2spin.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_spen/ufdosycosy_2spin.m`
- Signature: `ufdosycosy_2spin()`
- Total lines: 100

## Purpose

Ultrafast 3D DOSY-COSY for a two-spin system. Calculation time: hours on NVidia Tesla A100, much longer on CPU Ludmilla Guduff Jean-Nicolas Dumez

## Physical / mathematical content

- SPEN / ultrafast NMR examples. These files encode spatially dependent phase evolution and acquisition, linking pulse gradients, diffusion attenuation, and single-scan multidimensional encoding.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.
- Finite-difference discretisation appears in the implementation, so numerical accuracy depends on stencil order, boundary handling, and the balance between resolution and conditioning.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 10-11: Spin system; implemented by `sys.isotopes={'1H','1H'}`.
- Lines 13-14: Interactions; implemented by `sys.magnet=14.1`.
- Lines 19-20: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 23-24: Algorithmic options; implemented by `sys.disable={'pt'}`.
- Lines 27-28: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 31-32: Assumptions; implemented by `spin_system=assume(spin_system,'nmr')`.
- Lines 34-35: Sample geometry; implemented by `parameters.dims=0.015`.
- Lines 39-40: Relaxation phantom; implemented by `parameters.rlx_ph={}`.
- Lines 43-44: Initial and detection state phantoms; implemented by `parameters.rho0_ph={ones(parameters.npts,1)}`.
- Lines 49-50: Diffusion and flow; implemented by `parameters.u=zeros(parameters.npts,1)`.
- Lines 53-54: Acquisition parameters; implemented by `parameters.spins={'1H'}`.
- Lines 64-65: Third dimension parameters; implemented by `parameters.sweep=1302`.
- Lines 67-68: Encoding parameters; implemented by `parameters.pulsenpoints=1000`.
- Lines 76-77: Equation for the fit; implemented by `parameters.dscale=1e-10`.
- Lines 81-82: Window multiplication and apodization; implemented by `parameters.apodize='hamming'`.
- Lines 85-86: Fitting model; implemented by `parameters.model='keeler_corr'`.
- Lines 88-89: Simulation; implemented by `fid=imaging(spin_system,@spendosycosy,parameters)`.
- Lines 91-92: Processing; implemented by `spectrum=fftshift(fft(fid,[],1),1)`.

### Key state/data transformations

- Lines 11: computes `sys.isotopes` using `sys.isotopes={'1H','1H'}`.
- Lines 14: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 15: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={2.0,-1.3}`.
- Lines 16: computes `inter.coupling.scalar{1,2}` using `inter.coupling.scalar{1,2}=15`.
- Lines 17: computes `inter.coupling.scalar{2,2}` using `inter.coupling.scalar{2,2}=0`.
- Lines 20: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 21: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 24: computes `sys.disable` using `sys.disable={'pt'}`.
- Lines 25: computes `sys.enable` using `sys.enable={'greedy'}`.
- Lines 28: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 35: computes `parameters.dims` using `parameters.dims=0.015`.
- Lines 36: computes `parameters.npts` using `parameters.npts=3000`.
- Lines 37: computes `parameters.deriv` using `parameters.deriv={'period',7}`.
- Lines 40: computes `parameters.rlx_ph` using `parameters.rlx_ph={}`.
- Lines 41: computes `parameters.rlx_op` using `parameters.rlx_op={}`.
- Lines 44: computes `parameters.rho0_ph` using `parameters.rho0_ph={ones(parameters.npts,1)}`.
- Lines 45: computes `parameters.rho0_st` using `parameters.rho0_st={state(spin_system,'Lz','1H')}`.
- Lines 46: computes `parameters.coil_ph` using `parameters.coil_ph={ones(parameters.npts,1)}`.

## Implementation structure

- Ultrafast 3D DOSY-COSY for a two-spin system.
- Calculation time: hours on NVidia Tesla A100, much longer on CPU
- Ludmilla Guduff
- Jean-Nicolas Dumez
- Spin system
- Interactions
- Basis set
- Algorithmic options
- Spinach housekeeping
- Assumptions
- Sample geometry
- Relaxation phantom

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `assume()`, `state()`, `imaging()`, `fftshift()`, `kfigure()`, `volplot()`.
