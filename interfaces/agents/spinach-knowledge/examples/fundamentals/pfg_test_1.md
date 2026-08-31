# examples/fundamentals/pfg_test_1.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/fundamentals/pfg_test_1.m`
- Signature: `pfg_test_1()`
- Total lines: 77

## Purpose

A test of the explicit gradient pulse function that uses the auxiliary matrix formalism to compute sample volume integral. For details, see:

## Physical / mathematical content

- Fundamentals examples. These are unit tests, convention checks, and pedagogical demonstrations of operator algebra, perturbation theory, tensor conventions, symmetry, quadrature, and numerical differentiation.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 10-11: Magnet and isotopes; implemented by `sys.magnet=5.9`.
- Lines 14-15: Random chemical shifts and couplings; implemented by `inter.zeeman.scalar={10*rand(1),10*rand(1),10*rand(1)}`.
- Lines 21-22: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 25-26: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 30-31: Spin Hamiltonian; implemented by `L=hamiltonian(spin_system)`.
- Lines 33-34: Build initial state vector; implemented by `rho=rand(size(L,1),1)`.
- Lines 37-38: Determine projection quantum numbers of the basis; implemented by `[~,M]=lin2lm(spin_system.bas.basis)`.
- Lines 40-41: Determine the coherence order of each state; implemented by `coherence_orders=sum(M,2)`.
- Lines 43-44: Find out which coherence orders are present; implemented by `unique_coherence_orders=unique(coherence_orders).'`.
- Lines 46-47: Weight coherence orders by the number of states; implemented by `weighting=linspace(0.1,0.9,numel(unique_coherence_orders)); weight_idx=1`.
- Lines 54-55: Evolve under the homospoil pulse; implemented by `nsteps=100`.
- Lines 62-63: Preallocate the trajectory; implemented by `rho_stack=zeros([length(rho) nsteps],'like',1i)`.
- Lines 65-66: Compute the trajectory; implemented by `parfor t_idx=1:nsteps`.
- Lines 71-72: Analyze the trajectory; implemented by `kfigure()`.

### Control flow inferred from the code

- Line 48: `for` loop over `coherence_order=unique_coherence_orders`.
- Line 66: `parfor` loop over `t_idx=1:nsteps`.

### Key state/data transformations

- Lines 11: computes `sys.magnet` using `sys.magnet=5.9`.
- Lines 12: computes `sys.isotopes` using `sys.isotopes={'1H','1H','1H'}`.
- Lines 15: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={10*rand(1),10*rand(1),10*rand(1)}`.
- Lines 16: computes `inter.coupling.scalar` using `inter.coupling.scalar=cell(3)`.
- Lines 17: computes `inter.coupling.scalar{1,2}` using `inter.coupling.scalar{1,2}=20*rand(1)`.
- Lines 18: computes `inter.coupling.scalar{1,3}` using `inter.coupling.scalar{1,3}=20*rand(1)`.
- Lines 19: computes `inter.coupling.scalar{2,3}` using `inter.coupling.scalar{2,3}=20*rand(1)`.
- Lines 22: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 23: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 26: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 31: computes `L` using `L=hamiltonian(spin_system)`.
- Lines 34: computes `rho` using `rho=rand(size(L,1),1)`.
- Lines 35: computes `rho(1)` using `rho(1)=1; rho=rho./norm(rho)`.
- Lines 38: computes `[~,M]` using `[~,M]=lin2lm(spin_system.bas.basis)`.
- Lines 41: computes `coherence_orders` using `coherence_orders=sum(M,2)`.
- Lines 44: computes `unique_coherence_orders` using `unique_coherence_orders=unique(coherence_orders).'`.
- Lines 47: computes `weighting` using `weighting=linspace(0.1,0.9,numel(unique_coherence_orders)); weight_idx=1`.
- Lines 50: computes `rho(subspace_mask)` using `rho(subspace_mask)=rho(subspace_mask)./norm(rho(subspace_mask))*weighting(weight_idx)`.

## Implementation structure

- A test of the explicit gradient pulse function that uses the auxiliary
- matrix formalism to compute sample volume integral. For details, see:
- Magnet and isotopes
- Random chemical shifts and couplings
- Basis set
- Spinach housekeeping
- Spin Hamiltonian
- Build initial state vector
- Determine projection quantum numbers of the basis
- Determine the coherence order of each state
- Find out which coherence orders are present
- Weight coherence orders by the number of states

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `assume()`, `hamiltonian()`, `rho()`, `lin2lm()`, `weighting()`, `time_axis()`, `rho_stack()`, `grad_pulse()`, `kfigure()`, `trajan()`, `set()`.
