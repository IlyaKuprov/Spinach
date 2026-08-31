# examples/imaging/echo_planar_2d.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/imaging/echo_planar_2d.m`
- Signature: `echo_planar_2d()`
- Total lines: 90

## Purpose

Echo planar imaging example in 2D for a brain phantom. Simulation time: seconds, faster with a Tesla V100 GPU.

## Physical / mathematical content

- MRI and spectroscopic-imaging examples. These files combine gradient terms, spatial encoding, diffusion, slice selection, k-space sampling, and Fourier reconstruction, generally within Fokker-Planck or explicit spatial-grid descriptions.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 10-11: Isotopes; implemented by `sys.isotopes={'1H'}`.
- Lines 13-14: Magnetic induction; implemented by `sys.magnet=5.9`.
- Lines 16-17: Chemical shifts; implemented by `inter.zeeman.scalar={0.0}`.
- Lines 19-20: Relaxation model; implemented by `inter.relaxation={'t1_t2'}`.
- Lines 26-27: Disable path tracing; implemented by `sys.disable={'pt'}`.
- Lines 29-33: This needs a GPU sys.enable={'gpu'};; implemented by `bas.formalism='sphten-liouv'`.
- Lines 32-33: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 36-37: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 40-41: Sequence parameters; implemented by `parameters.spins={'1H'}`.
- Lines 50-51: Relaxation superoperators; implemented by `[R1Op,R2Op]=rlx_t1_t2(spin_system)`.
- Lines 53-54: Phantom library call; implemented by `[R1Ph,R2Ph,PDPh,dims,npts]=phantoms('brain-medres')`.
- Lines 57-58: Sample geometry; implemented by `parameters.dims=dims([1 2])`.
- Lines 62-63: Relaxation phantoms; implemented by `parameters.rlx_ph={R1Ph,R2Ph}`.
- Lines 66-67: Initial and detection state phantoms; implemented by `parameters.rho0_ph={PDPh}`.
- Lines 72-73: No diffusion or flow; implemented by `parameters.u=zeros(parameters.npts)`.
- Lines 77-78: Run the simulation; implemented by `mri=imaging(spin_system,@epi_2d,parameters)`.
- Lines 80-81: For FOV calculation, G{1} is effectively halved; implemented by `parameters.pe_grad_amp=parameters.pe_grad_amp/2`.
- Lines 83-84: Plotting; implemented by `kfigure(); scale_figure([2.5 1.0])`.

### Key state/data transformations

- Lines 11: computes `sys.isotopes` using `sys.isotopes={'1H'}`.
- Lines 14: computes `sys.magnet` using `sys.magnet=5.9`.
- Lines 17: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={0.0}`.
- Lines 20: computes `inter.relaxation` using `inter.relaxation={'t1_t2'}`.
- Lines 21: computes `inter.rlx_keep` using `inter.rlx_keep='diagonal'`.
- Lines 22: computes `inter.equilibrium` using `inter.equilibrium='zero'`.
- Lines 23: computes `inter.r1_rates` using `inter.r1_rates={1.0}`.
- Lines 24: computes `inter.r2_rates` using `inter.r2_rates={1.0}`.
- Lines 27: computes `sys.disable` using `sys.disable={'pt'}`.
- Lines 33: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 34: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 37: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 41: computes `parameters.spins` using `parameters.spins={'1H'}`.
- Lines 42: computes `parameters.decouple` using `parameters.decouple={}`.
- Lines 43: computes `parameters.offset` using `parameters.offset=0.0`.
- Lines 44: computes `parameters.image_size` using `parameters.image_size=[101 105]`.
- Lines 45: computes `parameters.ro_grad_amp` using `parameters.ro_grad_amp=5.3e-3`.
- Lines 46: computes `parameters.pe_grad_amp` using `parameters.pe_grad_amp=4.8e-3`.

## Implementation structure

- Echo planar imaging example in 2D for a brain phantom.
- Simulation time: seconds, faster with a Tesla V100 GPU.
- Isotopes
- Magnetic induction
- Chemical shifts
- Relaxation model
- Disable path tracing
- This needs a GPU
- sys.enable={'gpu'};
- Basis set
- Spinach housekeeping
- Sequence parameters

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `rlx_t1_t2()`, `phantoms()`, `R1Ph()`, `R2Ph()`, `PDPh()`, `dims()`, `npts()`, `state()`, `imaging()`, `kfigure()`, `scale_figure()`, `subplot()`, `mri_2d_plot()`, `ktitle()`.
