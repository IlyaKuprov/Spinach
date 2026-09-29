# examples/nmr_diffusion/flow_test_1b.m

[Source](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_diffusion/flow_test_1b.m)

This one-dimensional example transports a concentration profile by diffusion and flow without spin dynamics. Its source describes a periodic sample with a diffusion coefficient that increases quadratically from left to right; the flow is spatially constant. The source header estimates calculation time as seconds.

The sample dimension is `0.02`, with 100 points and derivative setting `{'period',7}`. The flow field is `u=0.3*ones(100,1)`. The spatial diffusion profile is `dxx=5e-5*(linspace(0,1,100).^2)'`, ranging from zero at the first grid point to `5e-5` at the last. No units are assigned in the source to these parameter values or to the timestep. The initial profile is `rho=exp(-0.125*((1:100)-20).^2)'`.

After constructing the ghost-spin system and basis, the example builds and inflates the transport generator with `v2fplanck` and propagates using `timestep=5e-4`, `nsteps=70`, and `evolution(spin_system,F,[],rho,timestep,nsteps,'trajectory')`. It plots each trajectory column against the sample coordinate in metres, with concentration labelled on the vertical axis and display limits `[0,1]`.