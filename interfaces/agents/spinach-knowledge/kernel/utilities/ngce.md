# kernel/utilities/ngce.m

**Source:** <https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/ngce.m>

## Purpose

Numerical integral route to the Redfield relaxation superoperator. The function computes a laboratory-frame relaxation superoperator `R` directly from a molecular dynamics trajectory of stochastic Hamiltonian superoperators, using numerical evaluation of Redfield's time integral, and optionally returns the element-by-element standard deviation of the mean of `R`.

Only a single chemical substance is supported. Segmented input raises
`Spinach:ngce:segmentedSubstances` before integration; the scalar unit-state
projection is not a direct-sum projector, including when regularisation is used.

## Numerical method and sampling regime

The zero-mean stochastic superoperators in `H1` are correlated across the molecular-dynamics trajectory while `H0` supplies the coherent propagator. A trapezium-rule lag integral over the estimated correlation time `tau_est` is averaged across trajectory stripes to form the real symmetric laboratory-frame Redfield relaxation superoperator. The unit-state component is protected from damping; optional `reg` regularises very small rates, and requesting `dR` returns an elementwise uncertainty of the mean across stripes. The result retains non-secular terms: any secular approximation is the caller’s responsibility.

This numerical estimate requires resolved dynamics rather than merely positive input values: at least 50 MD time steps per shortest `H0` period (`2*pi/normest(H0)`), at least 10 steps per `tau_est`, and a trajectory lasting at least 200 correlation-time estimates. A coarser step or shorter trajectory is rejected before integration.

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

The geometric unit-state projector is built with concentration one, independently of the initial population, including a zero-population single substance.
