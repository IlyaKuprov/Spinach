# examples/fundamentals/symmetry_4.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/fundamentals/symmetry_4.m`
- Signature: `symmetry_4()`
- Total lines: 50

## Purpose

Hamiltonian symmetrization for a radical pair with four equivalent nuclei under the S4 permutation group.

## Physical / mathematical content

- Fundamentals examples. These are unit tests, convention checks, and pedagogical demonstrations of operator algebra, perturbation theory, tensor conventions, symmetry, quadrature, and numerical differentiation.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 8-9: Magnetic field; implemented by `sys.magnet=0`.
- Lines 11-12: Spin system; implemented by `sys.isotopes ={'E','E','1H','1H','1H','1H'}`.
- Lines 14-15: Basis set; implemented by `bas.formalism='zeeman-hilb'`.
- Lines 20-21: Interactions; implemented by `inter.zeeman.scalar={2.002 2.002 0 0 0 0}`.
- Lines 28-29: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 32-33: Assumptions; implemented by `spin_system=assume(spin_system,'labframe')`.
- Lines 35-36: Hamiltonian superoperator; implemented by `H=hamiltonian(spin_system)`.
- Lines 38-39: Symmetry factorization; implemented by `S=horzcat(spin_system.bas.irrep.projector)`.
- Lines 41-42: Plotting; implemented by `kfigure(); scale_figure([1.5 1])`.

### Key state/data transformations

- Lines 9: computes `sys.magnet` using `sys.magnet=0`.
- Lines 12: computes `sys.isotopes` using `sys.isotopes ={'E','E','1H','1H','1H','1H'}`.
- Lines 15: computes `bas.formalism` using `bas.formalism='zeeman-hilb'`.
- Lines 16: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 17: computes `bas.sym_spins` using `bas.sym_spins={[3 4 5 6]}`.
- Lines 18: computes `bas.sym_group` using `bas.sym_group={'S4'}`.
- Lines 21: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={2.002 2.002 0 0 0 0}`.
- Lines 22: computes `inter.coupling.scalar` using `inter.coupling.scalar=num2cell(mt2hz([0 0 0.295 0.295 0.295 0.295`.
- Lines 29: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 36: computes `H` using `H=hamiltonian(spin_system)`.
- Lines 39: computes `S` using `S=horzcat(spin_system.bas.irrep.projector)`.

## Implementation structure

- Hamiltonian symmetrization for a radical pair with four
- equivalent nuclei under the S4 permutation group.
- Magnetic field
- Spin system
- Basis set
- Interactions
- Spinach housekeeping
- Assumptions
- Hamiltonian superoperator
- Symmetry factorization
- Plotting

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `num2cell()`, `mt2hz()`, `create()`, `basis()`, `assume()`, `hamiltonian()`, `horzcat()`, `kfigure()`, `scale_figure()`, `subplot()`, `spy()`, `ktitle()`, `xline()`, `yline()`.
