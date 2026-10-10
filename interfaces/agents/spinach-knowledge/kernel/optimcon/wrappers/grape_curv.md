# kernel/optimcon/wrappers/grape_curv.m

- Signature: `[traj_data,fidelity,df_du]=grape_curv(waveform_u,u2x,dx_du,spin_system)`

Evaluates the Cartesian GRAPE objective in user-defined coordinates. `waveform_u` has coordinates in rows and time samples in columns. At each time sample, `u2x` maps a coordinate column to control-operator coefficients; `dx_du` supplies the Jacobian with curvilinear coordinates along rows and Cartesian controls along columns, so the gradient pullback is `dx_du(u)*df_dx`.

With two requested outputs the wrapper returns trajectory and fidelity only. With three, it pulls every Cartesian gradient slice—including penalty slices—back to curvilinear coordinates and then zeros entries selected by `control.freeze`. The mask must match `waveform_u`. Penalties are evaluated in the rectilinear representation; this wrapper returns no Hessian.

Before calling `grape_xy`, the wrapper saves `control.freeze` locally and clears it so the Cartesian calculation does not apply a mask in the wrong coordinates. A nonempty mask must match `waveform_u` and cannot be combined with `control.basis`. Other guards require `spin_system.control` and function handles for `u2x` and `dx_du`, and `waveform_u` must be a real numeric array. Requests other than two or three outputs are rejected.

[Source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/optimcon/wrappers/grape_curv.m) · [Wiki](https://spindynamics.org/wiki/index.php?title=grape_curv.m)
