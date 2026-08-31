# examples/spin_chemistry/cidnp_geminate.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/spin_chemistry/cidnp_geminate.m`
- Signature: `cidnp_geminate()`
- Total lines: 59

## Purpose

A basic example of the geminate CIDNP effect simulation. Calculation time: seconds

## Physical / mathematical content

- Spin-chemistry examples. These scripts treat radical pairs, recombination channels, chemically induced dynamic nuclear polarisation, and magnetic-field effects. The theory combines spin-selective kinetics with singlet-triplet interconversion.
- The relevant state manifold is the singlet/triplet decomposition, where permutation symmetry controls selection rules, relaxation susceptibility, and convertibility to ordinary magnetisation.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 11-12: System specification; implemented by `sys.magnet=14.1`.
- Lines 21-22: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 25-26: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 29-30: Get the Hamiltonian; implemented by `H=hamiltonian(assume(spin_system,'esr'))`.
- Lines 32-33: Get the kinetics superoperator; implemented by `K=kinetics(spin_system)`.
- Lines 35-36: Get the initial state; implemented by `rho=singlet(spin_system,1,2)`.
- Lines 38-39: Double up the problem (no dynamics assumed in the product subspace); implemented by `H=[1*H 0*H; 0*H 0*H]`.
- Lines 42-43: Set up a reaction ("whatever is leaving reactants must appear in products"); implemented by `K=[1*K 0*K; -1*K 0*K]`.
- Lines 45-46: Assemble the Liouvillian; implemented by `L=H+1i*K`.
- Lines 48-49: Evolve for a microsecond; implemented by `rho=evolution(spin_system,L,[],rho,1e-6,1,'final')`.
- Lines 51-52: Check the nuclear magnetisation; implemented by `Nz=state(spin_system,'Lz','1H')`.

### Key state/data transformations

- Lines 12: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 13: computes `sys.isotopes` using `sys.isotopes ={'E','E','1H'}`.
- Lines 14: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={2.0023 2.0024 1.0}`.
- Lines 15: computes `inter.coupling.scalar` using `inter.coupling.scalar=cell(3,3)`.
- Lines 16: computes `inter.coupling.scalar{2,3}` using `inter.coupling.scalar{2,3}=1e7`.
- Lines 17: computes `inter.chem.rp_theory` using `inter.chem.rp_theory='haberkorn'`.
- Lines 18: computes `inter.chem.rp_electrons` using `inter.chem.rp_electrons=[1 2]`.
- Lines 19: computes `inter.chem.rp_rates` using `inter.chem.rp_rates=[1e7 0]`.
- Lines 22: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 23: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 26: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 30: computes `H` using `H=hamiltonian(assume(spin_system,'esr'))`.
- Lines 33: computes `K` using `K=kinetics(spin_system)`.
- Lines 36: computes `rho` using `rho=singlet(spin_system,1,2)`.
- Lines 46: computes `L` using `L=H+1i*K`.
- Lines 52: computes `Nz` using `Nz=state(spin_system,'Lz','1H')`.
- Lines 53: computes `rho_reac` using `rho_reac=rho(1:(numel(rho)/2))`.
- Lines 54: computes `rho_prod` using `rho_prod=rho((numel(rho)/2+1):end)`.

## Implementation structure

- A basic example of the geminate CIDNP effect simulation.
- Calculation time: seconds
- System specification
- Basis set
- Spinach housekeeping
- Get the Hamiltonian
- Get the kinetics superoperator
- Get the initial state
- Double up the problem (no dynamics assumed in the product subspace)
- Set up a reaction ("whatever is leaving reactants must appear in products")
- Assemble the Liouvillian
- Evolve for a microsecond

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `hamiltonian()`, `assume()`, `kinetics()`, `singlet()`, `evolution()`, `state()`, `rho()`, `num2str()`.
