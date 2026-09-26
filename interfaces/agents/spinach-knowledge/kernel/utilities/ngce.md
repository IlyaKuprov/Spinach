# kernel/utilities/ngce.m

- Signature: `[R,dR]=ngce(spin_system,H0,H1,dt,tau_est,reg)`

## Purpose

Numerical integral route to the Redfield relaxation superopera- tor. Syntax: [R,dR]=ngce(spin_system,H0,H1,dt,tau_est,reg)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- The relevant state manifold is the singlet/triplet decomposition, where permutation symmetry controls selection rules, relaxation susceptibility, and convertibility to ordinary magnetisation.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- The code contains an inverse-problem or ill-conditioning aspect and therefore introduces explicit regularisation, model selection, or stabilisation logic.

## Parameters / inputs

- H0 -static laboratory frame Hamiltonian commutation su-
- peroperator acting in the background, a matrix
- H1 -stochastic part (zero mean) of the laboratory frame
- Hamiltonian commutation superoperator, a cell array
- of matrices, one for each step of the MD trajectory.
- dt -time step of the MD trajectory, seconds
- tau_est -H1 autocorrelation time estimate for internal
- safety control, seconds
- reg -optional overall relaxation rate, this is added to
- every eigenvalue of the resulting matrix to prevent
- very small relaxation rates (e.g. singlets) from
- jumping into positive due to integration accuracy
- limits and then causing problems

## Outputs

- R -laboratory frame relaxation superoperator
- dR -standard deviation of the mean of R, element by element
- Note: enough trajectory points must be present to converge
- the ensemble averages and Redfield's integral.
- Note: the result is returned in the LABORATORY FRAME -eli-
- minating non-secular terms is user's responsibility.

## Implementation structure

- Numerical integral route to the Redfield relaxation superopera-
- tor. Syntax:
- [R,dR]=ngce(spin_system,H0,H1,dt,tau_est,reg)
- H0 -static laboratory frame Hamiltonian commutation su-
- peroperator acting in the background, a matrix
- H1 -stochastic part (zero mean) of the laboratory frame
- Hamiltonian commutation superoperator, a cell array
- of matrices, one for each step of the MD trajectory.
- dt -time step of the MD trajectory, seconds
- tau_est -H1 autocorrelation time estimate for internal
- safety control, seconds
- reg -optional overall relaxation rate, this is added to
