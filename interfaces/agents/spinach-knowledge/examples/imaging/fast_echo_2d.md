# examples/imaging/fast_echo_2d.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/imaging/fast_echo_2d.m`
- Signature: `fast_echo_2d()`
- Total lines: 83

## Purpose

Fast (in the experiment duration sense) spin echo 2D brain imaging example. Simulation time: hours, faster with a Tesla V100 GPU.

## Physical / mathematical content

- MRI and spectroscopic-imaging examples. These files combine gradient terms, spatial encoding, diffusion, slice selection, k-space sampling, and Fourier reconstruction, generally within Fokker-Planck or explicit spatial-grid descriptions.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 11-12: Isotopes; implemented by `sys.isotopes={'1H'}`.
- Lines 14-15: Magnetic induction; implemented by `sys.magnet=5.9`.
- Lines 17-18: Chemical shifts; implemented by `inter.zeeman.scalar={0.0}`.
- Lines 20-21: Relaxation model; implemented by `inter.relaxation={'t1_t2'}`.
- Lines 27-28: Disable path tracing; implemented by `sys.disable={'pt'}`.
- Lines 30-34: This needs a GPU sys.enable={'gpu'};; implemented by `bas.formalism='sphten-liouv'`.
- Lines 33-34: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 37-38: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 41-42: Sequence parameters; implemented by `parameters.spins={'1H'}`.
- Lines 51-52: Relaxation superoperators; implemented by `[R1Op,R2Op]=rlx_t1_t2(spin_system)`.
- Lines 54-55: Phantom library call; implemented by `[R1Ph,R2Ph,PDPh,dims,npts]=phantoms('brain-medres')`.
- Lines 58-59: Sample geometry; implemented by `parameters.dims=dims([1 2])`.
- Lines 63-64: Relaxation phantoms; implemented by `parameters.rlx_ph={R1Ph,R2Ph}`.
- Lines 67-68: Initial and detection state phantoms; implemented by `parameters.rho0_ph={PDPh}`.
- Lines 73-74: Run the simulation; implemented by `mri=imaging(spin_system,@fse,parameters)`.
- Lines 76-77: Plotting; implemented by `kfigure(); scale_figure([2.5 1.0])`.

### Key state/data transformations

- Lines 12: computes `sys.isotopes` using `sys.isotopes={'1H'}`.
- Lines 15: computes `sys.magnet` using `sys.magnet=5.9`.
- Lines 18: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={0.0}`.
- Lines 21: computes `inter.relaxation` using `inter.relaxation={'t1_t2'}`.
- Lines 22: computes `inter.rlx_keep` using `inter.rlx_keep='diagonal'`.
- Lines 23: computes `inter.equilibrium` using `inter.equilibrium='zero'`.
- Lines 24: computes `inter.r1_rates` using `inter.r1_rates={1}`.
- Lines 25: computes `inter.r2_rates` using `inter.r2_rates={1}`.
- Lines 28: computes `sys.disable` using `sys.disable={'pt'}`.
- Lines 34: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 35: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 38: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 42: computes `parameters.spins` using `parameters.spins={'1H'}`.
- Lines 43: computes `parameters.decouple` using `parameters.decouple={}`.
- Lines 44: computes `parameters.offset` using `parameters.offset=0.0`.
- Lines 45: computes `parameters.image_size` using `parameters.image_size=[101 105]`.
- Lines 46: computes `parameters.ro_grad_amp` using `parameters.ro_grad_amp=5.3e-3`.
- Lines 47: computes `parameters.ro_grad_dur` using `parameters.ro_grad_dur=2e-3`.

## Implementation structure

- Fast (in the experiment duration sense) spin echo 2D brain
- imaging example.
- Simulation time: hours, faster with a Tesla V100 GPU.
- Isotopes
- Magnetic induction
- Chemical shifts
- Relaxation model
- Disable path tracing
- This needs a GPU
- sys.enable={'gpu'};
- Basis set
- Spinach housekeeping

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `rlx_t1_t2()`, `phantoms()`, `R1Ph()`, `R2Ph()`, `PDPh()`, `dims()`, `npts()`, `state()`, `imaging()`, `kfigure()`, `scale_figure()`, `subplot()`, `mri_2d_plot()`, `ktitle()`.
