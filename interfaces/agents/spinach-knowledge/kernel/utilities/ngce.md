# kernel/utilities/ngce.m

- Signature: `[R,dR]=ngce(spin_system,H0,H1,dt,tau_est,reg)`

## Purpose

Computes a Redfield relaxation superoperator by numerical integration along a molecular-dynamics (MD) trajectory.

## Parameters

- `spin_system` — Spinach system object used for reporting, propagation, tolerances, and unit-state operators.
- `H0` — static laboratory-frame Hamiltonian commutation superoperator acting in the background; a matrix.
- `H1` — zero-mean stochastic part of the laboratory-frame Hamiltonian commutation superoperator; a cell array of matrices, one per MD trajectory step.
- `dt` — MD trajectory time step, in seconds.
- `tau_est` — estimate of the `H1` autocorrelation time, in seconds, used for internal safety checks.
- `reg` — optional overall relaxation rate. When nonzero, it is subtracted through `unit_oper(spin_system)` to prevent very small relaxation rates, such as singlet rates, from becoming positive because of integration accuracy limits.

## Outputs

- `R` — laboratory-frame relaxation superoperator.
- `dR` — element-by-element standard deviation of the mean of `R`, calculated when a second output is requested.

## Numerical method and checks

The routine estimates the shortest coherent-dynamics period as `2*pi/normest(H0)`. It requires at least 50 trajectory points per `H0` period, at least 10 points per `tau_est`, and a trajectory duration of at least `200*tau_est`; otherwise it raises an error. Enough trajectory points must be available to converge both the ensemble averages and Redfield's integral.

The correlation-time integration uses `5*ceil(tau_est/dt)` steps. Propagators are generated from `H0` and `dt` and cleaned using `spin_system.tols.prop_chop`. The trajectory is divided into complete, non-overlapping stripes of this length. A `parfor` loop applies the trapezium rule to each stripe's integral, with intermediate terms cleaned using `spin_system.tols.liouv_zero`. Each stripe contribution is restricted to its real symmetric, trace-preserving part before the contributions are averaged into `R`. When requested, `dR` is obtained from the element-by-element variance across stripes. After optional regularisation, the unit-state projector is used again so that the unit state is not damped.

`R` is returned in the **laboratory frame**. Eliminating non-secular terms is the user's responsibility.

## Reference and contacts

- <https://spindynamics.org/wiki/index.php?title=ngce.m>
- ilya.kuprov@weizmann.ac.il
- jpresteg@uga.edu
