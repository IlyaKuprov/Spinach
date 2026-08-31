# examples/nmr_diffusion/flow_test_1b.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_diffusion/flow_test_1b.m`
- Signature: `flow_test_1b()`
- Total lines: 62

## Purpose

A standard diffusion and flow equation solver with no spin dynamics present and periodic boundary condition. Diffusion coefficient increases quadratically from left to right. Calculation time: seconds.

## Physical / mathematical content

- Diffusion examples. The dominant mathematics is diffusion or advection-diffusion PDE propagation, sometimes with additional spin phase accumulation under gradients.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 12-13: Ghost spin; implemented by `sys.magnet=0`.
- Lines 16-17: No spin interactions; implemented by `inter.zeeman.matrix=cell(1)`.
- Lines 20-21: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 24-25: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 28-29: Sample geometry; implemented by `parameters.dims=0.02`.
- Lines 33-34: Diffusion and flow parameters; implemented by `parameters.u=0.3*ones(parameters.npts,1)`.
- Lines 37-38: Diffusion and flow generator; implemented by `F=v2fplanck(spin_system,parameters); F=inflate(F)`.
- Lines 40-41: Initial condition; implemented by `rho=exp(-0.125*((1:100)-20).^2)'`.
- Lines 43-44: Timing parameters; implemented by `timestep=5e-4; nsteps=70`.
- Lines 46-47: System trajectory; implemented by `traj=evolution(spin_system,F,[],rho,timestep,nsteps,'trajectory')`.
- Lines 49-50: Physically correct X axis; implemented by `x_axis=linspace(-parameters.dims/2,parameters.dims/2,parameters.npts)`.
- Lines 52-53: Plotting; implemented by `kfigure()`.

### Control flow inferred from the code

- Line 54: `for` loop over `n=1:nsteps`.

### Key state/data transformations

- Lines 13: computes `sys.magnet` using `sys.magnet=0`.
- Lines 14: computes `sys.isotopes` using `sys.isotopes={'G'}`.
- Lines 17: computes `inter.zeeman.matrix` using `inter.zeeman.matrix=cell(1)`.
- Lines 18: computes `inter.coupling.matrix` using `inter.coupling.matrix=cell(1)`.
- Lines 21: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 22: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 25: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 29: computes `parameters.dims` using `parameters.dims=0.02`.
- Lines 30: computes `parameters.npts` using `parameters.npts=100`.
- Lines 31: computes `parameters.deriv` using `parameters.deriv={'period',7}`.
- Lines 34: computes `parameters.u` using `parameters.u=0.3*ones(parameters.npts,1)`.
- Lines 35: computes `parameters.dxx` using `parameters.dxx=5e-5*(linspace(0,1,parameters.npts).^2)'`.
- Lines 38: computes `F` using `F=v2fplanck(spin_system,parameters); F=inflate(F)`.
- Lines 41: computes `rho` using `rho=exp(-0.125*((1:100)-20).^2)'`.
- Lines 44: computes `timestep` using `timestep=5e-4; nsteps=70`.
- Lines 47: computes `traj` using `traj=evolution(spin_system,F,[],rho,timestep,nsteps,'trajectory')`.
- Lines 50: computes `x_axis` using `x_axis=linspace(-parameters.dims/2,parameters.dims/2,parameters.npts)`.

## Implementation structure

- A standard diffusion and flow equation solver with no spin
- dynamics present and periodic boundary condition. Diffusion
- coefficient increases quadratically from left to right.
- Calculation time: seconds.
- Ghost spin
- No spin interactions
- Basis set
- Spinach housekeeping
- Sample geometry
- Diffusion and flow parameters
- Diffusion and flow generator
- Initial condition

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `v2fplanck()`, `inflate()`, `evolution()`, `kfigure()`, `traj()`, `kxlabel()`, `kylabel()`, `pause()`.
