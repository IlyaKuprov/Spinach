# examples/nmr_diffusion/diffusion_test_2c.m

- Signature: `diffusion_test_2c()`

## Purpose

Evolve a phantom under a two-dimensional diffusion and flow generator with no spin interactions. The sample geometry uses periodic derivatives. The source comments estimate a calculation time of minutes.

## Physical / mathematical content

- Loads `R1` from `phantom_a.mat` as the initial state.
- Configures a ghost spin (`'G'`) with zero magnetization and no Zeeman or coupling interactions.
- Uses a 108-by-90 grid over dimensions `[0.02 0.02]`, with `parameters.deriv={'period',7}`.
- Sets both flow fields, `u` and `v`, to zero. Sets each diffusion tensor field, `dxx`, `dxy`, `dyx`, and `dyy`, to a constant `5e-5` across the grid.

## Numerical / algorithmic content

- Creates the spin system with the `sphten-liouv` formalism and `none` basis approximation.
- Constructs the diffusion and flow generator with `v2fplanck(spin_system,parameters)` and applies `inflate` to it.
- Calls `evolution` in `'trajectory'` mode with initial state `R1(:)`, time step `5e-4`, and 200 steps.

## Implementation structure

- Loads the phantom and configures a ghost spin without spin interactions.
- Creates the spin system and basis, then specifies the sample geometry, zero flow fields, and diffusion tensor fields.
- Builds the diffusion and flow generator and computes the trajectory.
- Plots each of the 200 trajectory columns as a 108-by-90 image, calling `drawnow` and pausing for `0.025` seconds between frames.
