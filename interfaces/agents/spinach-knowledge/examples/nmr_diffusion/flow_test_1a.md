# examples/nmr_diffusion/flow_test_1a.m

[Source](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_diffusion/flow_test_1a.m)

This is a one-dimensional diffusion-and-flow transport example with no spin dynamics. It uses a constant flow field and constant diffusion coefficient on a periodic sample; the source labels the plotted sample coordinate in metres. The source header estimates calculation time as seconds.

The sample dimension is `0.02`, discretised at 100 points, and the derivative operator is configured as `{'period',7}`. The velocity field is `u=0.3*ones(100,1)` and the diffusion field is `dxx=5e-5*ones(100,1)`. The source does not annotate units for velocity, diffusion coefficient, or timestep. The initial concentration profile is `rho=exp(-0.125*((1:100)-20).^2)'`.

The ghost-spin system and basis are passed to `v2fplanck`, and the inflated generator is propagated with `timestep=5e-4`, `nsteps=70`, and `evolution(spin_system,F,[],rho,timestep,nsteps,'trajectory')`. The plot uses `x=linspace(-dims/2,dims/2,100)`, labels the axes “sample coordinate, m” and “concentration”, and plots every trajectory column with bounds `[-dims/2,dims/2]` and `[0,1]`.