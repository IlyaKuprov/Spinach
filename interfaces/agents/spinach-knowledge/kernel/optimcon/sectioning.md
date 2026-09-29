# kernel/optimcon/sectioning.m

Signature: `[alpha,fx_1,gfx_1,exitflag,data]=sectioning(cost_function,a,b,x_0,fx_0,gfx_0,dir,data,spin_system)`


It returns scalar step `alpha`, objective value `fx_1`, column gradient `gfx_1`, updated workspace `data`, and `exitflag`. This is the bracket-refinement part of the line search; it does not choose the search direction or form a Hessian. `a` and `b` are bracket endpoint structures containing `alpha`, `fx` and `gfx` fields. `x_0` and the gradients/direction are real column vectors; `fx_0` is a real scalar. Trial points are `x_0+alpha*dir`, evaluated with `objeval`.

At entry, if `spin_system.control.freeze` is nonempty, the routine multiplies both `dir` and the starting gradient `gfx_0` by the complement of `freeze(:)`. The trial gradient returned by `objeval` is not explicitly masked in this function. This routine does not read a phase-cycle mask.

Each trial is selected by cubic interpolation within a contracted interval whose endpoints use `spin_system.control.ls_tau2` and `spin_system.control.ls_tau3`. Failed monotonicity or sufficient-increase checks move the upper bracket; a trial that passes them is accepted only when the curvature test passes. Numerically unresolved interpolation or bracket collapse returns the lower endpoint with `exitflag=-2`; a Wolfe-accepted step returns `exitflag=0`.

The routine reads these line-search fields rather than assigning their defaults. The `optimcon` wrapper initialises `ls_tau2=0.1`, `ls_tau3=0.5`, `ls_c1=1e-2` and `ls_c2=0.9`. The acceptance tests are supplied by `alpha_conds`; `ls_tau1` is a bracketing expansion setting, not used in this sectioning routine.

[Source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/optimcon/sectioning.m)
[Spinach Wiki](https://spindynamics.org/wiki/index.php?title=sectioning.m)
