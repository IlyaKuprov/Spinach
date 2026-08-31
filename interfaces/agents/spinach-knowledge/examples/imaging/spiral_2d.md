# examples/imaging/spiral_2d.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/imaging/spiral_2d.m`
- Signature: `spiral_2d()`
- Total lines: 85

## Purpose

Spiral K-space imaging example in 2D. Calculation time: minutes. Ahmed Allami Ilya Kuprov

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
- Lines 26-27: Algorithmic options; implemented by `sys.disable={'pt'}`.
- Lines 29-30: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 33-34: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 37-38: Sequence parameters; implemented by `parameters.spins={'1H'}`.
- Lines 47-48: Sample geometry; implemented by `parameters.dims=[0.30 0.25]`.
- Lines 52-53: Relaxation phantom; implemented by `[R1,R2]=rlx_t1_t2(spin_system)`.
- Lines 58-59: Initial and detection state phantoms; implemented by `parameters.rho0_ph={ones(prod(parameters.npts,1))}`.
- Lines 64-65: Diffusion and flow; implemented by `parameters.u=zeros(parameters.npts)`.
- Lines 69-70: Run the simulation; implemented by `mri=imaging(spin_system,@spiral,parameters)`.
- Lines 72-73: Translate gradients for plotting; implemented by `parameters.ro_grad_dur=parameters.spiral_dur/sqrt(parameters.spiral_npts)/2`.
- Lines 78-79: Plotting; implemented by `loc=get(0,'defaultfigureposition'); figure('Position',[loc(1:2) 3*loc(3) loc(4)])`.

### Key state/data transformations

- Lines 11: computes `sys.isotopes` using `sys.isotopes={'1H'}`.
- Lines 14: computes `sys.magnet` using `sys.magnet=5.9`.
- Lines 17: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={0.0}`.
- Lines 20: computes `inter.relaxation` using `inter.relaxation={'t1_t2'}`.
- Lines 21: computes `inter.rlx_keep` using `inter.rlx_keep='diagonal'`.
- Lines 22: computes `inter.equilibrium` using `inter.equilibrium='zero'`.
- Lines 23: computes `inter.r1_rates` using `inter.r1_rates={30}`.
- Lines 24: computes `inter.r2_rates` using `inter.r2_rates={70}`.
- Lines 27: computes `sys.disable` using `sys.disable={'pt'}`.
- Lines 30: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 31: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 34: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 38: computes `parameters.spins` using `parameters.spins={'1H'}`.
- Lines 39: computes `parameters.decouple` using `parameters.decouple={}`.
- Lines 40: computes `parameters.offset` using `parameters.offset=0.0`.
- Lines 41: computes `parameters.grad_amp` using `parameters.grad_amp=5e-2`.
- Lines 42: computes `parameters.spiral_dur` using `parameters.spiral_dur=5e-3`.
- Lines 43: computes `parameters.spiral_npts` using `parameters.spiral_npts=3000`.

## Implementation structure

- Spiral K-space imaging example in 2D.
- Calculation time: minutes.
- Ahmed Allami
- Ilya Kuprov
- Isotopes
- Magnetic induction
- Chemical shifts
- Relaxation model
- Algorithmic options
- Basis set
- Spinach housekeeping
- Sequence parameters

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `rlx_t1_t2()`, `load()`, `state()`, `imaging()`, `get()`, `figure()`, `loc()`, `subplot()`, `mri_2d_plot()`.
