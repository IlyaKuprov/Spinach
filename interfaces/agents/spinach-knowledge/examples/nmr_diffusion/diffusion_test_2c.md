# examples/nmr_diffusion/diffusion_test_2c.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_diffusion/diffusion_test_2c.m`
- Signature: `diffusion_test_2c()`
- Total lines: 63

## Purpose

A standard diffusion equation solver with no spin dynamics present. Anisotropic diffusion with a pe- riodic boundary condition. Calculation time: minutes.

## Physical / mathematical content

- Diffusion examples. The dominant mathematics is diffusion or advection-diffusion PDE propagation, sometimes with additional spin phase accumulation under gradients.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 12-13: Load the phantom; implemented by `load('phantom_a.mat','R1')`.
- Lines 15-16: Ghost spin; implemented by `sys.magnet=0`.
- Lines 19-20: No spin interactions; implemented by `inter.zeeman.matrix=cell(1)`.
- Lines 23-24: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 27-28: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 31-32: Sample geometry; implemented by `parameters.dims=[0.02 0.02]`.
- Lines 36-37: 2D flow parameters; implemented by `parameters.u=zeros(parameters.npts)`.
- Lines 40-41: 2D diffusion tensor field; implemented by `parameters.dxx=5e-5*ones(parameters.npts)`.
- Lines 46-47: Diffusion and flow generator; implemented by `F=v2fplanck(spin_system,parameters); F=inflate(F)`.
- Lines 49-50: Timing parameters; implemented by `timestep=5e-4; nsteps=200`.
- Lines 52-53: System trajectory; implemented by `traj=evolution(spin_system,F,[],R1(:),timestep,nsteps,'trajectory')`.
- Lines 55-56: Plotting; implemented by `kfigure()`.

### Control flow inferred from the code

- Line 57: `for` loop over `n=1:nsteps`.

### Key state/data transformations

- Lines 16: computes `sys.magnet` using `sys.magnet=0`.
- Lines 17: computes `sys.isotopes` using `sys.isotopes={'G'}`.
- Lines 20: computes `inter.zeeman.matrix` using `inter.zeeman.matrix=cell(1)`.
- Lines 21: computes `inter.coupling.matrix` using `inter.coupling.matrix=cell(1)`.
- Lines 24: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 25: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 28: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 32: computes `parameters.dims` using `parameters.dims=[0.02 0.02]`.
- Lines 33: computes `parameters.npts` using `parameters.npts=[108 90]`.
- Lines 34: computes `parameters.deriv` using `parameters.deriv={'period',7}`.
- Lines 37: computes `parameters.u` using `parameters.u=zeros(parameters.npts)`.
- Lines 38: computes `parameters.v` using `parameters.v=zeros(parameters.npts)`.
- Lines 41: computes `parameters.dxx` using `parameters.dxx=5e-5*ones(parameters.npts)`.
- Lines 42: computes `parameters.dxy` using `parameters.dxy=5e-5*ones(parameters.npts)`.
- Lines 43: computes `parameters.dyx` using `parameters.dyx=5e-5*ones(parameters.npts)`.
- Lines 44: computes `parameters.dyy` using `parameters.dyy=5e-5*ones(parameters.npts)`.
- Lines 47: computes `F` using `F=v2fplanck(spin_system,parameters); F=inflate(F)`.
- Lines 50: computes `timestep` using `timestep=5e-4; nsteps=200`.

## Implementation structure

- A standard diffusion equation solver with no spin
- dynamics present. Anisotropic diffusion with a pe-
- riodic boundary condition.
- Calculation time: minutes.
- Load the phantom
- Ghost spin
- No spin interactions
- Basis set
- Spinach housekeeping
- Sample geometry
- 2D flow parameters
- 2D diffusion tensor field

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `load()`, `create()`, `basis()`, `v2fplanck()`, `inflate()`, `evolution()`, `kfigure()`, `traj()`, `pause()`.
