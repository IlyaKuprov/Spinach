# kernel/optimcon/tgrape.m

Signature: `[fidelity,grad]=tgrape(spin_system,drift,controls,waveform,dt_grid,time_unit,rho_init,rho_targ)`


This GRAPE objective propagates an initial Liouville-space column state through the control slices and returns `real(rho_targ'*fwd_traj(:,end))`. The accepted formalisms are `sphten-liouv` and `zeeman-liouv`. `drift` and every cell in `controls` must be same-size square matrices; `rho_init` and `rho_targ` are compatible column vectors. `waveform` is a real numeric array with one row per control and one column per slice; its coefficients are in rad/s.

`dt_grid` is a finite real column vector with one entry per waveform column. The source checks finiteness but does not require each duration entry to be positive. `time_unit` is a finite positive real scalar in seconds. Before propagation, the function replaces `dt_grid` with `dt_grid*time_unit`; after computing each slice derivative it multiplies the returned gradient by `time_unit`. `grad` is a column vector with one entry per slice, giving derivatives with respect to the original `dt_grid` coordinates. The source documentation recommends scaling the duration variables so entries of `dt_grid` are of order 1; this is a scaling rationale, not an enforced bound.

For slice `n`, the gradient is the real part of `-1i*bwd_traj(:,n)'*L*fwd_traj(:,n)`, where `L=drift+sum_k waveform(k,n)*controls{k}`. No freeze or phase-cycle mask is read by this routine, and it supplies no Hessian.

[Source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/optimcon/tgrape.m)
[Spinach Wiki](https://spindynamics.org/wiki/index.php?title=tgrape.m)
