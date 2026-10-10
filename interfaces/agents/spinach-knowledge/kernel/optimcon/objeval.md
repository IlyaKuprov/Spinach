# kernel/optimcon/objeval.m

- Signature: `[data,fx,grad,hess]=objeval(x,objfun_handle,data,spin_system)`

## Purpose

Adapts an objective-function call to the number of outputs requested by an optimisation routine, and combines the objective's fidelity, gradient, and Hessian components. The source notes that this function will be eliminated in a future release.

## Inputs and guards

All four arguments are required; the source assigns no defaults. It checks that `objfun_handle` is a function handle and that `x` is non-empty, numeric, and real. The guard does not explicitly test `isvector(x)`, despite the error text describing a vector. The routine also expects `data.x_shape` and the counter fields it updates; it does not validate `data` or `spin_system` here.

## Wrapper transformations

The wrapper reshapes `x` to `data.x_shape` and calls the objective handle with that shaped control and `spin_system`. With two outputs requested from `objeval`, it asks the objective for `[traj_data,fidelity]`; with three, it also requests a gradient; with four, it requests a Hessian as well. Other output counts raise an error.

In all supported cases, the combined objective is `fidelity(1)-sum(fidelity(2:end))`: the first fidelity component is added and the remaining components are subtracted. For gradient calls, the same first-component-minus-rest combination is applied along the third dimension, then the result is flattened to a column with `grad(:)`. For Hessian calls the component combination is applied along the third dimension and the resulting matrix is returned without flattening. The routine stores the uncombined values in `data.fx_sep_pen` and the trajectory data in `data.traj_data`.

Each supported call increments `data.count.fx`; gradient requests also increment `data.count.gfx`, and Hessian requests increment `data.count.hfx`. This wrapper does not select line-search or Hessian-update options and does not apply freeze or phase-cycle masks.

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/optimcon/objeval.m)
[Spinach Wiki page](https://spindynamics.org/wiki/index.php?title=objeval.m)
