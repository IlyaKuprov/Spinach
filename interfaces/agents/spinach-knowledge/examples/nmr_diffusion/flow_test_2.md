# examples/nmr_diffusion/flow_test_2.m

- Signature: `flow_test_2()`

## Purpose

Demonstrates combined diffusion and flow in two dimensions with periodic boundary conditions. The source notes a calculation time of minutes.

## Physical / mathematical content

- A loaded phantom, `R1`, evolves under a two-dimensional flow and diffusion model. The flow components are uniform (`u=v=0.2`); the diffusion tensor has diagonal components `dxx=dyy=5e-5` and zero cross-components.
- The sample dimensions are `[0.02 0.02]`, with a `[108 90]` grid and periodic derivatives specified by `{'period',7}`. The ghost spin has no spin interactions.

## Numerical / algorithmic content

- `v2fplanck(spin_system,parameters)` constructs the diffusion-and-flow generator, which is then passed through `inflate`.
- `evolution` computes the trajectory from `R1(:)` using a timestep of `5e-4` for 200 steps. The trajectory frames are reshaped to `[108 90]` and displayed.

## Implementation structure

- Load `R1` from `phantom_a.mat`.
- Configure a ghost spin with no spin interactions and create the Spinach system using the `sphten-liouv` formalism with no basis approximation.
- Set the sample geometry, periodic derivatives, two-dimensional flow field, and two-dimensional diffusion tensor field.
- Build the diffusion-and-flow generator with `v2fplanck`, compute the loaded phantom’s trajectory with `evolution`, and display its frames.
