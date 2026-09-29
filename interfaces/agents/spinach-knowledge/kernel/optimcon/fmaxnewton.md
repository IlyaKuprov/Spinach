# kernel/optimcon/fmaxnewton.m

[Source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/optimcon/fmaxnewton.m) · [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=fmaxnewton.m)

## Purpose and interface

`[x,data]=fmaxnewton(spin_system,cost_function,guess)` seeks a local maximum using the method selected in `spin_system.control.method`. The prepared `spin_system` comes from `optimcon.m`; the objective is a function handle, and `guess` is real numeric. The search vector is flattened internally and the returned `x` is reshaped to `size(guess)`. `data.x_shape` records that shape. The counters `data.count.iter`, `fx`, `gfx`, `hfx`, and `rfo` track iterations, objective, gradient, Hessian, and RFO evaluations; an objective may add other `data` fields.

## Search methods and derivatives

The supported methods are `lbfgs`, `rbfgs`, `newton`, and `goodwin`. On the first iteration, the objective and gradient are evaluated; `newton` and `goodwin` also request a Hessian each iteration. Their Hessian is symmetrised and made real, restricted to unfrozen coordinates, then regularised by `hessreg` before the direction is computed. `goodwin` uses the regularised Newton method with Goodwin acceleration.

`lbfgs` forms a limited-memory direction from recent step and gradient differences without Hessian regularisation. `rbfgs` forms a BFGS pseudo-Hessian from the same kind of history and regularises it before solving for a direction. Histories are capped by `spin_system.control.n_grads`. Both start from a scaled gradient direction. The optimiser rejects non-finite directions and falls back to the projected gradient if a direction is not an ascent direction.

## Line search and stopping

Each iteration brackets an acceptable step with `bracketing`; when that routine requests it, `sectioning` searches the resulting interval. There is no separate line-search selector in this function. A failed sectioning search (exit flag `-2`) prevents the step update. Otherwise the step is accepted, and a configured checkpoint is saved.

The first-iteration guard rejects an assembled gradient whose 2-norm on unfrozen coordinates is below `1e-6`. A step 1-norm below `tol_x` or an unfrozen-gradient 2-norm below `tol_g` terminates the search. `optimcon.m` supplies defaults of `max_iter=100`, `tol_x=1e-3`, `tol_g=1e-6`, and `n_grads=50` (the history length is used by `lbfgs` and `rbfgs`). With `max_iter=0`, the routine evaluates only the objective; it performs no iteration and does not apply the initial-gradient guard.

## Masks, waveform shape, and checkpointing

An empty `freeze` mask leaves all coordinates free. Otherwise its dimensions must exactly match `guess`; the mask is flattened and frozen coordinates are excluded from the search direction and gradient tests. `optimcon.m` rejects a mask that freezes the entire waveform. A configured checkpoint is loaded from the global scratch directory when the file exists, replacing the supplied initial point after a `numel` consistency check.

For the `rectangle` integrator, `guess` has `pulse_nsteps` columns, or one column per basis coefficient when a basis is used; the basis has `pulse_nsteps` columns. For `trapezium`, those widths are `pulse_nsteps+1`; a basis must have that many columns. Other integrators are rejected. When `video_file` is configured, the routine writes Motion JPEG AVI frames during iterations.
