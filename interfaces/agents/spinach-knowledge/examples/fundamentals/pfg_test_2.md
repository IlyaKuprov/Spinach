# examples/fundamentals/pfg_test_2.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/fundamentals/pfg_test_2.m`
- Signature: `pfg_test_2()`
- Total lines: 86

## Purpose

Demonstrate the use of the auxiliary matrix algorithm in generating a gradient sandwich multiple-quantum filter. For further details see: Calculation time: seconds

## Physical / mathematical content

- Fundamentals examples. These are unit tests, convention checks, and pedagogical demonstrations of operator algebra, perturbation theory, tensor conventions, symmetry, quadrature, and numerical differentiation.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 12-13: Magnet and isotopes; implemented by `sys.magnet=5.9`.
- Lines 16-17: Random chemical shifts and couplings; implemented by `inter.zeeman.scalar={10*rand(1),10*rand(1),10*rand(1)}`.
- Lines 23-24: Set the basis; implemented by `bas.formalism='sphten-liouv'`.
- Lines 27-28: Run Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 32-33: Build the Hamiltonian; implemented by `H=hamiltonian(spin_system)`.
- Lines 35-36: Propagator for a pi/2 pulse; implemented by `Hx=operator(spin_system,'Lx','1H')`.
- Lines 39-40: Build initial state vector; implemented by `rho=rand(size(H,1),1)`.
- Lines 43-44: Determine projection quantum numbers of the basis; implemented by `[~,M]=lin2lm(spin_system.bas.basis)`.
- Lines 46-47: Determine the coherence order of each state; implemented by `coherence_orders=sum(M,2)`.
- Lines 49-50: Find out which coherence orders are present; implemented by `unique_coherence_orders=unique(coherence_orders).'`.
- Lines 52-53: Choose initial state weightings; implemented by `weighting=[0,0,0,0,0,1,0]`.
- Lines 55-56: Weight coherence orders by the number of states; implemented by `weight_idx=1`.
- Lines 63-64: Build coherence diagram; implemented by `gradient_strengths=20*[1 1]`.
- Lines 70-71: Preallocate the trajectory; implemented by `rho_stack=zeros([length(rho) nsteps],'like',1i)`.
- Lines 73-74: Compute the trajectory; implemented by `parfor time_idx=1:nsteps`.
- Lines 80-81: Analyze the trajectory; implemented by `kfigure()`.

### Control flow inferred from the code

- Line 57: `for` loop over `coherence_order=unique_coherence_orders`.
- Line 74: `parfor` loop over `time_idx=1:nsteps`.

### Key state/data transformations

- Lines 13: computes `sys.magnet` using `sys.magnet=5.9`.
- Lines 14: computes `sys.isotopes` using `sys.isotopes={'1H','1H','1H'}`.
- Lines 17: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={10*rand(1),10*rand(1),10*rand(1)}`.
- Lines 18: computes `inter.coupling.scalar` using `inter.coupling.scalar=cell(3)`.
- Lines 19: computes `inter.coupling.scalar{1,2}` using `inter.coupling.scalar{1,2}=20*rand(1)`.
- Lines 20: computes `inter.coupling.scalar{1,3}` using `inter.coupling.scalar{1,3}=20*rand(1)`.
- Lines 21: computes `inter.coupling.scalar{2,3}` using `inter.coupling.scalar{2,3}=20*rand(1)`.
- Lines 24: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 25: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 28: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 33: computes `H` using `H=hamiltonian(spin_system)`.
- Lines 36: computes `Hx` using `Hx=operator(spin_system,'Lx','1H')`.
- Lines 37: computes `P` using `P=propagator(spin_system,Hx,pi/2)`.
- Lines 40: computes `rho` using `rho=rand(size(H,1),1)`.
- Lines 41: computes `rho(1)` using `rho(1)=1; rho=rho./norm(rho)`.
- Lines 44: computes `[~,M]` using `[~,M]=lin2lm(spin_system.bas.basis)`.
- Lines 47: computes `coherence_orders` using `coherence_orders=sum(M,2)`.
- Lines 50: computes `unique_coherence_orders` using `unique_coherence_orders=unique(coherence_orders).'`.

## Implementation structure

- Demonstrate the use of the auxiliary matrix algorithm in generating a
- gradient sandwich multiple-quantum filter. For further details see:
- Calculation time: seconds
- Magnet and isotopes
- Random chemical shifts and couplings
- Set the basis
- Run Spinach housekeeping
- Build the Hamiltonian
- Propagator for a pi/2 pulse
- Build initial state vector
- Determine projection quantum numbers of the basis
- Determine the coherence order of each state

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `assume()`, `hamiltonian()`, `operator()`, `propagator()`, `rho()`, `lin2lm()`, `weighting()`, `durations_temp()`, `time_axis()`, `rho_stack()`, `grad_sandw()`, `kfigure()`, `trajan()`, `set()`.
