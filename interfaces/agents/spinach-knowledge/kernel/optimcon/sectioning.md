# kernel/optimcon/sectioning.m

- Signature: `[alpha,fx_1,gfx_1,exitflag,data]=sectioning(cost_function,a,b,x_0,fx_0,...`

```matlab
[alpha,fx_1,gfx_1,exitflag,data]=sectioning(cost_function,a,b,x_0,fx_0,gfx_0,dir,data,spin_system)
```

Refines a step bracket by cubic interpolation until Wolfe conditions pass or numerical resolution is lost. `cost_function` is an objective handle; `a,b` are bracket structs with `alpha,fx,gfx`; `x_0,fx_0,gfx_0` are the current point, value and gradient; `dir` is the search direction; `data` is the workspace; `spin_system` supplies line-search settings. Frozen coordinates are masked. Trials are evaluated by `objeval`; Armijo, monotonicity and curvature tests update or accept the bracket. Returns step, value, gradient, updated workspace and `exitflag` (0 success, -2 failure).

[Source](https://spindynamics.org/wiki/index.php?title=sectioning.m)