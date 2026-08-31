# examples/nmr_spen/ufdosy_1spin.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_spen/ufdosy_1spin.m`
- Signature: `ufdosy_1spin()`
- Total lines: 87

## Purpose

Ultrafast DOSY for one spin. Calculation time: seconds on NVidia Tesla A100, much longer on CPU Ludmilla Guduff Jean-Nicolas Dumez

## Physical / mathematical content

- SPEN / ultrafast NMR examples. These files encode spatially dependent phase evolution and acquisition, linking pulse gradients, diffusion attenuation, and single-scan multidimensional encoding.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.
- Finite-difference discretisation appears in the implementation, so numerical accuracy depends on stencil order, boundary handling, and the balance between resolution and conditioning.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 10-11: Spin system; implemented by `sys.isotopes={'1H'}`.
- Lines 13-14: Interactions; implemented by `sys.magnet=14.1`.
- Lines 17-18: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 21-22: Algorithmic options; implemented by `sys.disable={'pt'}`.
- Lines 25-26: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 29-30: Assumptions; implemented by `spin_system=assume(spin_system,'nmr')`.
- Lines 32-33: Sample geometry; implemented by `parameters.dims=0.015`.
- Lines 37-38: Relaxation phantom; implemented by `parameters.rlx_ph={}`.
- Lines 41-42: Initial and detection state phantoms; implemented by `parameters.rho0_ph={ones(parameters.npts,1)}`.
- Lines 47-48: Diffusion and flow; implemented by `parameters.u=zeros(parameters.npts,1)`.
- Lines 51-52: Acquisition parameters; implemented by `parameters.spins={'1H'}`.
- Lines 60-61: Encoding parameters; implemented by `parameters.pulsenpoints=1000`.
- Lines 69-70: Simulation; implemented by `fid=imaging(spin_system,@spendosy,parameters)`.
- Lines 72-73: Processing; implemented by `squarespectrum=fftshift(fft(fid,[],1),1)`.
- Lines 76-77: Plotting; implemented by `npoints2=size(squarespectrum,2)`.

### Key state/data transformations

- Lines 11: computes `sys.isotopes` using `sys.isotopes={'1H'}`.
- Lines 14: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 15: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={7.0}`.
- Lines 18: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 19: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 22: computes `sys.disable` using `sys.disable={'pt'}`.
- Lines 23: computes `sys.enable` using `sys.enable={'greedy'}`.
- Lines 26: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 33: computes `parameters.dims` using `parameters.dims=0.015`.
- Lines 34: computes `parameters.npts` using `parameters.npts=3000`.
- Lines 35: computes `parameters.deriv` using `parameters.deriv={'period',7}`.
- Lines 38: computes `parameters.rlx_ph` using `parameters.rlx_ph={}`.
- Lines 39: computes `parameters.rlx_op` using `parameters.rlx_op={}`.
- Lines 42: computes `parameters.rho0_ph` using `parameters.rho0_ph={ones(parameters.npts,1)}`.
- Lines 43: computes `parameters.rho0_st` using `parameters.rho0_st={state(spin_system,'Lz','1H')}`.
- Lines 44: computes `parameters.coil_ph` using `parameters.coil_ph={ones(parameters.npts,1)}`.
- Lines 45: computes `parameters.coil_st` using `parameters.coil_st={state(spin_system,'L+','1H')}`.
- Lines 48: computes `parameters.u` using `parameters.u=zeros(parameters.npts,1)`.

## Implementation structure

- Ultrafast DOSY for one spin.
- Calculation time: seconds on NVidia Tesla A100, much longer on CPU
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

- Called routines detected from the main body: `create()`, `basis()`, `assume()`, `state()`, `imaging()`, `fftshift()`, `spin()`, `kfigure()`, `kxlabel()`, `kylabel()`.
