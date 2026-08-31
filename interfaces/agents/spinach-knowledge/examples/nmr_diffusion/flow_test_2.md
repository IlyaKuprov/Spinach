# examples/nmr_diffusion/flow_test_2.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_diffusion/flow_test_2.m`
- Signature: `flow_test_2()`
- Total lines: 62

## Purpose

A combination of diffusion and flow in two dimensions with a periodic boundary condition. Calculation time: minutes.

## Physical / mathematical content

- Diffusion examples. The dominant mathematics is diffusion or advection-diffusion PDE propagation, sometimes with additional spin phase accumulation under gradients.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 11-12: Load the phantom; implemented by `load('phantom_a.mat','R1')`.
- Lines 14-15: Ghost spin; implemented by `sys.magnet=0`.
- Lines 18-19: No spin interactions; implemented by `inter.zeeman.matrix=cell(1)`.
- Lines 22-23: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 26-27: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 30-31: Sample geometry; implemented by `parameters.dims=[0.02 0.02]`.
- Lines 35-36: 2D flow field; implemented by `parameters.u=0.2*ones(parameters.npts)`.
- Lines 39-40: 2D diffusion tensor field; implemented by `parameters.dxx=5e-5*ones(parameters.npts)`.
- Lines 45-46: Diffusion and flow generator; implemented by `F=v2fplanck(spin_system,parameters); F=inflate(F)`.
- Lines 48-49: Timing parameters; implemented by `timestep=5e-4; nsteps=200`.
- Lines 51-52: System trajectory; implemented by `traj=evolution(spin_system,F,[],R1(:),timestep,nsteps,'trajectory')`.
- Lines 54-55: Plotting; implemented by `kfigure()`.

### Control flow inferred from the code

- Line 56: `for` loop over `n=1:nsteps`.

### Key state/data transformations

- Lines 15: computes `sys.magnet` using `sys.magnet=0`.
- Lines 16: computes `sys.isotopes` using `sys.isotopes={'G'}`.
- Lines 19: computes `inter.zeeman.matrix` using `inter.zeeman.matrix=cell(1)`.
- Lines 20: computes `inter.coupling.matrix` using `inter.coupling.matrix=cell(1)`.
- Lines 23: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 24: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 27: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 31: computes `parameters.dims` using `parameters.dims=[0.02 0.02]`.
- Lines 32: computes `parameters.npts` using `parameters.npts=[108 90]`.
- Lines 33: computes `parameters.deriv` using `parameters.deriv={'period',7}`.
- Lines 36: computes `parameters.u` using `parameters.u=0.2*ones(parameters.npts)`.
- Lines 37: computes `parameters.v` using `parameters.v=0.2*ones(parameters.npts)`.
- Lines 40: computes `parameters.dxx` using `parameters.dxx=5e-5*ones(parameters.npts)`.
- Lines 41: computes `parameters.dxy` using `parameters.dxy=zeros(parameters.npts)`.
- Lines 42: computes `parameters.dyx` using `parameters.dyx=zeros(parameters.npts)`.
- Lines 43: computes `parameters.dyy` using `parameters.dyy=5e-5*ones(parameters.npts)`.
- Lines 46: computes `F` using `F=v2fplanck(spin_system,parameters); F=inflate(F)`.
- Lines 49: computes `timestep` using `timestep=5e-4; nsteps=200`.

## Implementation structure

- A combination of diffusion and flow in two dimensions
- with a periodic boundary condition.
- Calculation time: minutes.
- Load the phantom
- Ghost spin
- No spin interactions
- Basis set
- Spinach housekeeping
- Sample geometry
- 2D flow field
- 2D diffusion tensor field
- Diffusion and flow generator

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `load()`, `create()`, `basis()`, `v2fplanck()`, `inflate()`, `evolution()`, `kfigure()`, `traj()`, `pause()`.
