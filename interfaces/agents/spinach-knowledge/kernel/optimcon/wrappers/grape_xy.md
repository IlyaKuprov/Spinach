# kernel/optimcon/wrappers/grape_xy.m

- Signature: `[traj_data,fidelity,grad,hess]=grape_xy(waveform,spin_system)`

Evaluates the ensemble GRAPE objective for Cartesian control amplitudes and adds configured penalty terms using their weights and lower/upper bounds. This routine computes objective values and requested derivatives for a supplied waveform; it does not run the optimiser or choose a line-search option.

`waveform` must be a real numeric array with `ncontrols` rows. With an empty `control.basis`, its columns match `pulse_ntpts`. With a basis, the basis has `pulse_ntpts` columns and waveform has one column per basis function; the physical waveform is `waveform*control.basis`. A nonempty `control.freeze` cannot be combined with a basis; otherwise the wrapper passes `spin_system` through to `ensemble`.

Two outputs return trajectory and a fidelity row containing the ensemble objective followed by one weighted slice per penalty. Three outputs add `grad`, with one third-dimension slice per objective or penalty. Four outputs also add `hess`, whose first two dimensions are `numel(waveform)` and whose third dimension has the same slice convention. When a basis is used, gradients are pulled back to basis coefficients and Hessians are transformed on both coordinate dimensions. There is no explicit default optimiser, line-search, or Hessian-approximation setting in this wrapper.

[Source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/optimcon/wrappers/grape_xy.m) · [Wiki](https://spindynamics.org/wiki/index.php?title=grape_xy.m)
