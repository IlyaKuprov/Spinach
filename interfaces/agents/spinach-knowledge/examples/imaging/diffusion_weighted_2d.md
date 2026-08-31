# examples/imaging/diffusion_weighted_2d.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/imaging/diffusion_weighted_2d.m`
- Signature: `diffusion_weighted_2d()`
- Total lines: 80

## Purpose

2D diffusion weighted image with an arbitrary geometric pattern serving as diffusion coefficient distribution. Simulation time: minutes, faster with a Tesla V100 GPU.

## Physical / mathematical content

- MRI and spectroscopic-imaging examples. These files combine gradient terms, spatial encoding, diffusion, slice selection, k-space sampling, and Fourier reconstruction, generally within Fokker-Planck or explicit spatial-grid descriptions.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 11-12: Isotopes; implemented by `sys.isotopes={'1H'}`.
- Lines 14-15: Magnetic induction; implemented by `sys.magnet=5.9`.
- Lines 17-18: Chemical shifts; implemented by `inter.zeeman.scalar={0.0}`.
- Lines 20-24: This needs a GPU sys.enable={'gpu'};; implemented by `bas.formalism='sphten-liouv'`.
- Lines 23-24: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 27-28: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 31-32: Sequence parameters; implemented by `parameters.spins={'1H'}`.
- Lines 43-44: Sample geometry; implemented by `parameters.dims=[0.30 0.25]`.
- Lines 48-49: Relaxation phantoms and operators; implemented by `parameters.rlx_ph={}; parameters.rlx_op={}`.
- Lines 51-52: Initial and detection state phantoms; implemented by `parameters.rho0_ph={ones(prod(parameters.npts,1))}`.
- Lines 57-58: Diffusion and flow; implemented by `parameters.u=zeros(parameters.npts)`.
- Lines 61-62: 2D diffusion tensor field; implemented by `load('../../etc/phantoms/pattern.mat','pattern')`.
- Lines 68-69: Run the simulation; implemented by `mri=imaging(spin_system,@phase_enc_2d,parameters)`.
- Lines 71-72: Plotting; implemented by `loc=get(0,'defaultfigureposition')`.

### Key state/data transformations

- Lines 12: computes `sys.isotopes` using `sys.isotopes={'1H'}`.
- Lines 15: computes `sys.magnet` using `sys.magnet=5.9`.
- Lines 18: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={0.0}`.
- Lines 24: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 25: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 28: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 32: computes `parameters.spins` using `parameters.spins={'1H'}`.
- Lines 33: computes `parameters.decouple` using `parameters.decouple={}`.
- Lines 34: computes `parameters.offset` using `parameters.offset=0.0`.
- Lines 35: computes `parameters.image_size` using `parameters.image_size=[101 105]`.
- Lines 36: computes `parameters.diff_g_amp` using `parameters.diff_g_amp=[1e-3 1e-3]`.
- Lines 37: computes `parameters.ro_grad_amp` using `parameters.ro_grad_amp=4.3e-3`.
- Lines 38: computes `parameters.ro_grad_dur` using `parameters.ro_grad_dur=2e-3`.
- Lines 39: computes `parameters.pe_grad_amp` using `parameters.pe_grad_amp=3.8e-3`.
- Lines 40: computes `parameters.pe_grad_dur` using `parameters.pe_grad_dur=1e-3`.
- Lines 41: computes `parameters.t_echo` using `parameters.t_echo=1e-2`.
- Lines 44: computes `parameters.dims` using `parameters.dims=[0.30 0.25]`.
- Lines 45: computes `parameters.npts` using `parameters.npts=[90 108]`.

## Implementation structure

- 2D diffusion weighted image with an arbitrary geometric
- pattern serving as diffusion coefficient distribution.
- Simulation time: minutes, faster with a Tesla V100 GPU.
- Isotopes
- Magnetic induction
- Chemical shifts
- This needs a GPU
- sys.enable={'gpu'};
- Basis set
- Spinach housekeeping
- Sequence parameters
- Sample geometry

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `load()`, `imaging()`, `get()`, `figure()`, `loc()`, `subplot()`, `mri_2d_plot()`, `ktitle()`.
