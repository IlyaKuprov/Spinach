# examples/nmr_diffusion/diffusion_test_1.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_diffusion/diffusion_test_1.m`
- Signature: `diffusion_test_1()`
- Total lines: 61

## Purpose

A standard diffusion equation solver with no spin dynamics present. Calculation time: seconds

## Physical / mathematical content

- Diffusion examples. The dominant mathematics is diffusion or advection-diffusion PDE propagation, sometimes with additional spin phase accumulation under gradients.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 11-12: Ghost spin; implemented by `sys.magnet=0`.
- Lines 15-16: No spin interactions; implemented by `inter.zeeman.matrix=cell(1)`.
- Lines 19-20: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 23-24: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 27-28: Sample geometry; implemented by `parameters.dims=0.02`.
- Lines 32-33: Diffusion and flow parameters; implemented by `parameters.u=zeros(parameters.npts,1)`.
- Lines 36-37: Diffusion and flow generator; implemented by `F=v2fplanck(spin_system,parameters); F=inflate(F)`.
- Lines 39-40: Initial condition; implemented by `rho=exp(-0.125*((1:100)-20).^2)'`.
- Lines 42-43: Timing parameters; implemented by `timestep=5e-4; nsteps=90`.
- Lines 45-46: System trajectory; implemented by `traj=evolution(spin_system,F,[],rho,timestep,nsteps,'trajectory')`.
- Lines 48-49: Physically correct X axis; implemented by `x_axis=linspace(-parameters.dims/2,parameters.dims/2,parameters.npts)`.
- Lines 51-52: Plotting; implemented by `kfigure()`.

### Control flow inferred from the code

- Line 53: `for` loop over `n=1:nsteps`.

### Key state/data transformations

- Lines 12: computes `sys.magnet` using `sys.magnet=0`.
- Lines 13: computes `sys.isotopes` using `sys.isotopes={'G'}`.
- Lines 16: computes `inter.zeeman.matrix` using `inter.zeeman.matrix=cell(1)`.
- Lines 17: computes `inter.coupling.matrix` using `inter.coupling.matrix=cell(1)`.
- Lines 20: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 21: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 24: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 28: computes `parameters.dims` using `parameters.dims=0.02`.
- Lines 29: computes `parameters.npts` using `parameters.npts=100`.
- Lines 30: computes `parameters.deriv` using `parameters.deriv={'period',7}`.
- Lines 33: computes `parameters.u` using `parameters.u=zeros(parameters.npts,1)`.
- Lines 34: computes `parameters.diff` using `parameters.diff=5e-5`.
- Lines 37: computes `F` using `F=v2fplanck(spin_system,parameters); F=inflate(F)`.
- Lines 40: computes `rho` using `rho=exp(-0.125*((1:100)-20).^2)'`.
- Lines 43: computes `timestep` using `timestep=5e-4; nsteps=90`.
- Lines 46: computes `traj` using `traj=evolution(spin_system,F,[],rho,timestep,nsteps,'trajectory')`.
- Lines 49: computes `x_axis` using `x_axis=linspace(-parameters.dims/2,parameters.dims/2,parameters.npts)`.

## Implementation structure

- A standard diffusion equation solver with no spin
- dynamics present.
- Calculation time: seconds
- Ghost spin
- No spin interactions
- Basis set
- Spinach housekeeping
- Sample geometry
- Diffusion and flow parameters
- Diffusion and flow generator
- Initial condition
- Timing parameters

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `v2fplanck()`, `inflate()`, `evolution()`, `kfigure()`, `traj()`, `kylabel()`, `kxlabel()`, `pause()`.
