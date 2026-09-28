# examples/nmr_diffusion/diffusion_test_1.m

- Signature: `diffusion_test_1()`

## Purpose

Demonstrates diffusion of an initial concentration profile with no spin dynamics. The source notes a calculation time of seconds.

## Physical / mathematical content

- A ghost spin is defined with no Zeeman or coupling interactions.
- The sample spans 0.02 m and is represented by 100 spatial points. The flow velocity is zero at every point, and the diffusion parameter is `5e-5`.
- The initial concentration is a Gaussian-shaped profile centered near spatial grid point 20. The calculated trajectory follows its evolution over 90 time steps of `5e-4` each.
- The plot’s x-axis is the physical sample coordinate in metres, running from `-0.01` to `0.01`; it is not a time axis.

## Numerical / algorithmic content

- The spatial derivative setting is `{'period',7}`. `v2fplanck` constructs the diffusion-and-flow generator, which is then passed through `inflate`.
- `evolution(spin_system,F,[],rho,timestep,nsteps,'trajectory')` returns the trajectory. The source does not specify how `evolution` implements propagation.

## Implementation structure

- A standard diffusion equation solver with no spin
- dynamics present.
- Calculation time: seconds
- Ghost spin: sets zero magnetic field and isotope `G`.
- No spin interactions: supplies empty Zeeman and coupling matrices.
- Basis set: selects `sphten-liouv` formalism with no approximation.
- Spinach housekeeping: creates the spin system and constructs its basis.
- Sample geometry: sets the extent, point count, and derivative setting.
- Diffusion and flow parameters: sets zero flow and the diffusion parameter.
- Diffusion and flow generator: constructs and inflates `F`.
- Initial condition: constructs the Gaussian-shaped column vector `rho`.
- Timing parameters: sets the time step and number of steps.
- System trajectory: calls `evolution` in `trajectory` mode.
- Physically correct X axis: uses 100 evenly spaced coordinates across the sample.
- Plotting: plots concentration against sample coordinate for each trajectory column, with axes fixed to the sample extent and concentrations from 0 to 1; each frame is drawn with a `0.025`-second pause.
