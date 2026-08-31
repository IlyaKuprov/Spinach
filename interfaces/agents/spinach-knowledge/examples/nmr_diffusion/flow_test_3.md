# examples/nmr_diffusion/flow_test_3.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_diffusion/flow_test_3.m`
- Signature: `flow_test_3()`
- Total lines: 81

## Purpose

Circular flow in three-dimensional space in the absence of spin dynamics. Calculation time: minutes, faster on GPU.

## Physical / mathematical content

- Diffusion examples. The dominant mathematics is diffusion or advection-diffusion PDE propagation, sometimes with additional spin phase accumulation under gradients.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 11-12: Ghost spin; implemented by `sys.isotopes={'G'}`.
- Lines 15-16: No spin interactions; implemented by `inter.zeeman.matrix=cell(1)`.
- Lines 19-20: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 23-24: Algorithmic options; implemented by `sys.enable={'greedy'}`.
- Lines 26-27: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 30-31: Sample geometry; implemented by `parameters.dims=[0.02 0.02 0.02]`.
- Lines 35-38: Get a 3D grid; implemented by `[X,Y,Z]=ndgrid(linspace(-parameters.dims(1)/2,parameters.dims(1)/2,parameters.npts(1)), linspace(-parameters.dims(2)/2,parameters.dims(2)/2,parameters.npts(2)), linspace…`.
- Lines 40-41: Get circular wind vectors; implemented by `parameters.u=-1000*Y`.
- Lines 45-46: Constant diffusion tensor field; implemented by `parameters.dxx=8e-6*ones(parameters.npts)`.
- Lines 56-57: Get the initial condition; implemented by `X0=[-0.003 0.001 0.003]`.
- Lines 65-66: Fokker-Planck evolution generator; implemented by `F=v2fplanck(spin_system,parameters); F=inflate(F)`.
- Lines 68-69: System trajectory; implemented by `traj=evolution(spin_system,F,[],A(:),5e-5,200,'trajectory'); traj=full(traj)`.
- Lines 71-72: Plotting; implemented by `kfigure()`.

### Control flow inferred from the code

- Line 73: `for` loop over `n=1:size(traj,2)`.

### Key state/data transformations

- Lines 12: computes `sys.isotopes` using `sys.isotopes={'G'}`.
- Lines 13: computes `sys.magnet` using `sys.magnet=5.9`.
- Lines 16: computes `inter.zeeman.matrix` using `inter.zeeman.matrix=cell(1)`.
- Lines 17: computes `inter.coupling.matrix` using `inter.coupling.matrix=cell(1)`.
- Lines 20: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 21: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 24: computes `sys.enable` using `sys.enable={'greedy'}`.
- Lines 27: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 31: computes `parameters.dims` using `parameters.dims=[0.02 0.02 0.02]`.
- Lines 32: computes `parameters.npts` using `parameters.npts=[50 50 50]`.
- Lines 33: computes `parameters.deriv` using `parameters.deriv={'period',7}`.
- Lines 36-38: computes `[X,Y,Z]` using `[X,Y,Z]=ndgrid(linspace(-parameters.dims(1)/2,parameters.dims(1)/2,parameters.npts(1)), linspace(-parameters.dims(2)/2,parameters.dims(2)/2,parameters.npts(2)), linspace…`.
- Lines 41: computes `parameters.u` using `parameters.u=-1000*Y`.
- Lines 42: computes `parameters.v` using `parameters.v=+1000*X`.
- Lines 43: computes `parameters.w` using `parameters.w=zeros(size(X))`.
- Lines 46: computes `parameters.dxx` using `parameters.dxx=8e-6*ones(parameters.npts)`.
- Lines 47: computes `parameters.dxy` using `parameters.dxy=zeros(parameters.npts)`.
- Lines 48: computes `parameters.dxz` using `parameters.dxz=zeros(parameters.npts)`.

## Implementation structure

- Circular flow in three-dimensional space in the absence
- of spin dynamics.
- Calculation time: minutes, faster on GPU.
- Ghost spin
- No spin interactions
- Basis set
- Algorithmic options
- Spinach housekeeping
- Sample geometry
- Get a 3D grid
- Get circular wind vectors
- Constant diffusion tensor field

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `v2fplanck()`, `inflate()`, `evolution()`, `kfigure()`, `volplot()`, `traj()`.
