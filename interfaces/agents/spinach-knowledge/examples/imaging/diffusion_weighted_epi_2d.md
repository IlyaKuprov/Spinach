# examples/imaging/diffusion_weighted_epi_2d.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/imaging/diffusion_weighted_epi_2d.m`
- Signature: `diffusion_weighted_epi_2d()`
- Total lines: 95

## Purpose

2D echo planar imaging example in the presence of istropic diffusion. Stejskal-Tanner SE echo planar diffusion-weighted pulse sequence from Figure 1 in (https://doi.org/10.1148/radiol.09090021). Simulation time: minutes, faster with a Tesla V100 GPU.

## Physical / mathematical content

- MRI and spectroscopic-imaging examples. These files combine gradient terms, spatial encoding, diffusion, slice selection, k-space sampling, and Fourier reconstruction, generally within Fokker-Planck or explicit spatial-grid descriptions.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 12-13: Isotopes; implemented by `sys.isotopes={'1H'}`.
- Lines 15-16: Magnetic induction; implemented by `sys.magnet=5.9`.
- Lines 18-19: Chemical shifts; implemented by `inter.zeeman.scalar={0.0}`.
- Lines 21-22: Relaxation model; implemented by `inter.relaxation={'t1_t2'}`.
- Lines 28-29: Disable path tracing; implemented by `sys.disable={'pt'}`.
- Lines 31-35: This needs a GPU sys.enable={'gpu'};; implemented by `bas.formalism='sphten-liouv'`.
- Lines 34-35: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 38-39: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 42-43: Sequence parameters; implemented by `parameters.spins={'1H'}`.
- Lines 54-55: Relaxation superoperators; implemented by `[R1Op,R2Op]=rlx_t1_t2(spin_system)`.
- Lines 57-58: Phantom library call; implemented by `[R1Ph,R2Ph,PDPh,dims,npts]=phantoms('brain-medres')`.
- Lines 61-62: Sample geometry; implemented by `parameters.dims=dims([1 2])`.
- Lines 66-67: Relaxation phantoms; implemented by `parameters.rlx_ph={R1Ph,R2Ph}`.
- Lines 70-71: Initial and detection state phantoms; implemented by `parameters.rho0_ph={PDPh}`.
- Lines 76-77: 2D diffusion tensor field; implemented by `parameters.dxx=1e-4*ones(parameters.npts)`.
- Lines 82-83: Run the simulation; implemented by `mri=imaging(spin_system,@epi_2d,parameters)`.
- Lines 85-86: For FOV calculation, G{1} is effectively halved; implemented by `parameters.pe_grad_amp=parameters.pe_grad_amp/2`.
- Lines 88-89: Plotting; implemented by `kfigure(); scale_figure([2.5 1.0])`.

### Key state/data transformations

- Lines 13: computes `sys.isotopes` using `sys.isotopes={'1H'}`.
- Lines 16: computes `sys.magnet` using `sys.magnet=5.9`.
- Lines 19: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={0.0}`.
- Lines 22: computes `inter.relaxation` using `inter.relaxation={'t1_t2'}`.
- Lines 23: computes `inter.rlx_keep` using `inter.rlx_keep='diagonal'`.
- Lines 24: computes `inter.equilibrium` using `inter.equilibrium='zero'`.
- Lines 25: computes `inter.r1_rates` using `inter.r1_rates={1.0}`.
- Lines 26: computes `inter.r2_rates` using `inter.r2_rates={1.0}`.
- Lines 29: computes `sys.disable` using `sys.disable={'pt'}`.
- Lines 35: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 36: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 39: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 43: computes `parameters.spins` using `parameters.spins={'1H'}`.
- Lines 44: computes `parameters.decouple` using `parameters.decouple={}`.
- Lines 45: computes `parameters.offset` using `parameters.offset=0.0`.
- Lines 46: computes `parameters.image_size` using `parameters.image_size=[101 105]`.
- Lines 47: computes `parameters.ro_grad_amp` using `parameters.ro_grad_amp=5.3e-3`.
- Lines 48: computes `parameters.pe_grad_amp` using `parameters.pe_grad_amp=4.8e-3`.

## Implementation structure

- 2D echo planar imaging example in the presence of istropic diffusion.
- Stejskal-Tanner SE echo planar diffusion-weighted pulse sequence from
- Figure 1 in (https://doi.org/10.1148/radiol.09090021).
- Simulation time: minutes, faster with a Tesla V100 GPU.
- Isotopes
- Magnetic induction
- Chemical shifts
- Relaxation model
- Disable path tracing
- This needs a GPU
- sys.enable={'gpu'};
- Basis set

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `rlx_t1_t2()`, `phantoms()`, `R1Ph()`, `R2Ph()`, `PDPh()`, `dims()`, `npts()`, `state()`, `imaging()`, `kfigure()`, `scale_figure()`, `subplot()`, `mri_2d_plot()`, `ktitle()`.
