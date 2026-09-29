# kernel/utilities/ngce.m

**Source:** <https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/ngce.m>

## Purpose

Numerical integral route to the Redfield relaxation superoperator. The function computes a laboratory-frame relaxation superoperator `R` directly from a molecular dynamics trajectory of stochastic Hamiltonian superoperators, using numerical evaluation of Redfield's time integral, and optionally returns the element-by-element standard deviation of the mean of `R`.

## Behaviour

- Syntax: `[R,dR]=ngce(spin_system,H0,H1,dt,tau_est,reg)`.
- Input consistency is enforced first by an internal `grumble` check: `H0` must be a matrix, `H1` must be a cell array of matrices, and both `dt` and `tau_est` must be positive real numbers.
- The coherent dynamics timescale is computed as `2*pi/normest(H0)`, i.e. the shortest period in `H0`.
- The number of trajectory points per correlation time estimate is `ceil(tau_est/dt)`, and the number of points under the tau integral is five times that value (`n_tau_int_steps=5*npts_in_tau_c`).
- Timing diagnostics are printed via `report`, including the shortest period in `H0`, trajectory step length, total trajectory length, total trajectory points, the user estimate for `tau_c`, trajectory points per `H0` period, trajectory points per `tau_est`, and trajectory duration divided by `tau_est`.
- Sampling requirements are enforced with hard errors:
  - `timescale/dt<50` raises `insufficient H0 dynamics sampling, reduce trajectory time step.`
  - `tau_est/dt<10` raises `insufficient correlation function sampling, reduce trajectory time step.`
  - `traj_dur/tau_est<200` raises `insufficient ensemble average sampling, increase trajectory duration.`
- Propagators for the tau integrals are computed with `propagator(spin_system,H0,dt)`; the cumulative products `P{n}=P_dt*P{n-1}` are sparsified with `clean_up` using `spin_system.tols.prop_chop`, with `P{1}` set to the sparse identity.
- The trajectory is cut into `floor(traj_npts/n_tau_int_steps)` stripes (ensemble instances), truncated to a whole number of stripes, and reshaped so each stripe contains `n_tau_int_steps` points.
- For each stripe, Redfield's integral is evaluated with the trapezium rule over tau, with integrand terms of the form `-dt*H1s{1}*Ps{tau}*H1s{tau}*Ps{tau}'` cleaned up using `spin_system.tols.liouv_zero`.
- Each stripe integral is symmetrised as `real(F+F')/2` and made trace-preserving by subtracting `(U'*F*U)*USP`, where `U=unit_state(spin_system)` and `USP=U*U'`.
- Stripe integrals are accumulated in parallel (`parfor`) as sparse triplets and assembled into `R_sum`; the superoperator is `R=R_sum/nstripes`, keeping the real symmetric part.
- Error analysis is performed only when a second output argument is requested (`nargout>1`): the variance of the mean is `(R_sq_sum-(R_sum.^2)/nstripes)/(nstripes*(nstripes-1))`, negative variances are clipped to zero, and `dR=sqrt(dR_var)` element by element.
- If the optional `reg` argument exists and is nonzero, regularisation subtracts `reg*unit_oper(spin_system)` from `R`, adding `reg` to every eigenvalue to prevent very small relaxation rates (e.g. singlets) from jumping into positive values due to integration accuracy limits.
- Finally, the unit state is not damped: `R=R-(U'*R*U)*USP`.
- The result is returned in the laboratory frame; eliminating non-secular terms is the user's responsibility.
- Enough trajectory points must be present to converge the ensemble averages and Redfield's integral.

## Inputs and outputs

**Inputs**

- `spin_system` — Spinach spin system object supplying tolerances (`spin_system.tols.prop_chop`, `spin_system.tols.liouv_zero`) and reporting.
- `H0` — static laboratory-frame Hamiltonian commutation superoperator acting in the background, a matrix.
- `H1` — stochastic part (zero mean) of the laboratory-frame Hamiltonian commutation superoperator, a cell array of matrices, one for each step of the MD trajectory.
- `dt` — time step of the MD trajectory, seconds.
- `tau_est` — `H1` autocorrelation time estimate for internal safety control, seconds.
- `reg` — optional overall relaxation rate added to every eigenvalue of the resulting matrix to prevent very small relaxation rates (e.g. singlets) from jumping into positive due to integration accuracy limits and then causing problems.

**Outputs**

- `R` — laboratory-frame relaxation superoperator.
- `dR` — standard deviation of the mean of `R`, element by element.

## References

1. Spinach Wiki page for `ngce.m`: <https://spindynamics.org/wiki/index.php?title=ngce.m>
2. Spinach GitHub source file: <https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/ngce.m>
