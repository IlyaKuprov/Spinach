# examples/nmr_diffusion/diffusion_test_2c.m

[Source](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_diffusion/diffusion_test_2c.m)

This example propagates a loaded 2D phantom through a diffusion-and-flow transport generator without spin interactions. It loads `R1` from `phantom_a.mat` as the initial state and sets the velocity fields `u` and `v` to zero. The source comment labels the case anisotropic diffusion with periodic boundary conditions; the configured tensor components are all equal. The source header estimates calculation time as minutes.

The grid has dimensions `[0.02 0.02]` and `[108 90]` points. The derivative setting is `{'period',7}`. Each diffusion-tensor component, `dxx`, `dxy`, `dyx`, and `dyy`, is set uniformly to `5e-5`; the source gives no units for these values or the grid dimensions. With a ghost spin and empty Zeeman and coupling matrices, the example constructs the system and basis, then builds the generator using `v2fplanck(spin_system,parameters)` and `inflate`.

It requests 200 trajectory steps at `timestep=5e-4` with `evolution(spin_system,F,[],R1(:),timestep,nsteps,'trajectory')`. Each column is reshaped to `108-by-90` and displayed with `imagesc`; the loop calls `drawnow` and pauses `0.025` seconds per frame.