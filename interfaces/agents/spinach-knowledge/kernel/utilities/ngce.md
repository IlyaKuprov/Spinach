# kernel/utilities/ngce.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/ngce.m`
- Signature: `[R,dR]=ngce(spin_system,H0,H1,dt,tau_est,reg)`
- Total lines: 193

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
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 45-46: Check consistency; implemented by `grumble(H0,H1,dt,tau_est)`.
- Lines 48-49: Coherent dynamics timescale; implemented by `timescale=2*pi/normest(H0)`.
- Lines 51-52: Number of points in tau_c; implemented by `npts_in_tau_c=ceil(tau_est/dt)`.
- Lines 54-55: Number of points under tau integral; implemented by `n_tau_int_steps=5*npts_in_tau_c`.
- Lines 57-58: Trajectory point count; implemented by `traj_npts=numel(H1)`.
- Lines 60-61: Trajectory duration; implemented by `traj_dur=dt*(traj_npts-1)`.
- Lines 63-64: Print timing diagnostics; implemented by `report(spin_system,' ')`.
- Lines 71-72: Enforce sufficient H0 sampling; implemented by `report(spin_system,['Trajectory points per H0 period: ' num2str(timescale/dt)])`.
- Lines 77-78: Enforce sufficient tau_c sampling; implemented by `report(spin_system,['Trajectory points per tau_est: ' num2str(tau_est/dt)])`.
- Lines 83-84: Enforce sufficient trajectory duration; implemented by `report(spin_system,['Trajectory duration / tau_est: ' num2str(dt*(numel(H1)-1)/tau_est)])`.
- Lines 89-90: Get propagators for tau integrals; implemented by `report(spin_system,'computing H0 exponentials ')`.
- Lines 97-98: Determine the number of stripes; implemented by `nstripes=floor(traj_npts/n_tau_int_steps)`.
- Lines 102-103: Cut the trajectory up into stripes; implemented by `H1=H1(1:nstripes*n_tau_int_steps)`.
- Lines 106-107: Get the unit state projector; implemented by `U=unit_state(spin_system); USP=U*U'`.
- Lines 109-110: Error analysis switch; implemented by `calc_err=nargout>1`.
- Lines 112-113: Compute trajectory stripe integrals; implemented by `report(spin_system,'computing Redfield''s integral ')`.
- Lines 118-119: Get the stripe started; implemented by `H1s=H1(:,s); Ps=P; F=sparse(0)`.
- Lines 121-123: Trapezium rule for tau integral; implemented by `f_curr=clean_up(spin_system,-dt*H1s{1}*Ps{1}*H1s{1}*Ps{1}', spin_system.tols.liouv_zero)`.

### Control flow inferred from the code

- Line 73: conditional branch on `timescale/dt<50`.
- Line 79: conditional branch on `tau_est/dt<10`.
- Line 85: conditional branch on `traj_dur/tau_est<200`.
- Line 93: `for` loop over `n=2:n_tau_int_steps`.
- Line 116: `parfor` loop over `s=1:nstripes`.
- Line 124: `for` loop over `tau=2:n_tau_int_steps`.
- Line 138: conditional branch on `calc_err, vals_sq{s}=stripe_vals.^2; end`.
- Line 148: conditional branch on `calc_err`.
- Line 160: conditional branch on `exist('reg','var')&&(reg~=0)`.

### Key state/data transformations

- Lines 49: computes `timescale` using `timescale=2*pi/normest(H0)`.
- Lines 52: computes `npts_in_tau_c` using `npts_in_tau_c=ceil(tau_est/dt)`.
- Lines 55: computes `n_tau_int_steps` using `n_tau_int_steps=5*npts_in_tau_c`.
- Lines 58: computes `traj_npts` using `traj_npts=numel(H1)`.
- Lines 61: computes `traj_dur` using `traj_dur=dt*(traj_npts-1)`.
- Lines 91: computes `P` using `P=cell(n_tau_int_steps,1); P{1}=speye(size(H0))`.
- Lines 92: computes `P_dt` using `P_dt=propagator(spin_system,H0,dt)`.
- Lines 94: computes `P{n}` using `P{n}=clean_up(spin_system,P_dt*P{n-1},spin_system.tols.prop_chop)`.
- Lines 98: computes `nstripes` using `nstripes=floor(traj_npts/n_tau_int_steps)`.
- Lines 103: computes `H1` using `H1=H1(1:nstripes*n_tau_int_steps)`.
- Lines 107: computes `U` using `U=unit_state(spin_system); USP=U*U'`.
- Lines 110: computes `calc_err` using `calc_err=nargout>1`.
- Lines 114: computes `rows` using `rows=cell(nstripes,1); cols=cell(nstripes,1)`.
- Lines 115: computes `vals` using `vals=cell(nstripes,1); vals_sq=cell(nstripes,1)`.
- Lines 119: computes `H1s` using `H1s=H1(:,s); Ps=P; F=sparse(0)`.
- Lines 122-123: computes `f_curr` using `f_curr=clean_up(spin_system,-dt*H1s{1}*Ps{1}*H1s{1}*Ps{1}', spin_system.tols.liouv_zero)`.
- Lines 125-126: computes `f_next` using `f_next=clean_up(spin_system,-dt*H1s{1}*Ps{tau}*H1s{tau}*Ps{tau}', spin_system.tols.liouv_zero)`.
- Lines 127: computes `F` using `F=F+(f_curr+f_next)/2; f_curr=f_next`.

### Local helper functions

- Line 170: `grumble()` — `function grumble(H0,H1,dt,tau_est)`.
  - Representative operation: `if ~isnumeric(H0)`.
  - Representative operation: `error('H0 must be a matrix.')`.

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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `normest()`, `report()`, `num2str()`, `speye()`, `propagator()`, `clean_up()`, `unit_state()`, `cell2mat()`, `dR_var()`, `exist()`, `unit_oper()`, `iscell()`.
