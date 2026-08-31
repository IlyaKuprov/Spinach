# examples/esr_sol_swept/fieldsweep_porphyrin.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/esr_sol_swept/fieldsweep_porphyrin.m`
- Signature: `fieldsweep_porphyrin()`
- Total lines: 69

## Purpose

Field swept EPR spectrum of copper porphyrin complex, computed by finding resonance fields and transition moments. Calculation time: minutes.

## Physical / mathematical content

- Field-swept ESR examples. These files emphasise resonance-field finding, powder averaging, anisotropic g and hyperfine tensors, and intensity accumulation over orientation manifolds.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 10-11: Magnet field; implemented by `sys.magnet=1`.
- Lines 13-14: Isotopes; implemented by `sys.isotopes={'14N','14N','14N','14N','E','63Cu'}`.
- Lines 16-17: Array preallocation; implemented by `inter.zeeman.eigs=cell(6,1)`.
- Lines 23-24: Zeeman interactions; implemented by `inter.zeeman.eigs{5,1}=[2.0509 2.0509 2.1801]`.
- Lines 27-28: Hyperfine interactions; implemented by `inter.coupling.eigs{5,6}=[-70.9257 -70.9257 -575.0219]*1e6`.
- Lines 35-36: Basis set; implemented by `bas.formalism='zeeman-hilb'`.
- Lines 39-40: Symmetry; implemented by `bas.sym_group={'S4'}`.
- Lines 43-44: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 47-48: Experiment parameters; implemented by `parameters.spins={'E'}`.
- Lines 58-59: Run the simulation in the high-T approximation; implemented by `parameters.rho0=-state(spin_system,'Lz','E')`.
- Lines 62-63: Plotting; implemented by `kfigure(); plot(parameters.b_axis,spec)`.

### Key state/data transformations

- Lines 11: computes `sys.magnet` using `sys.magnet=1`.
- Lines 14: computes `sys.isotopes` using `sys.isotopes={'14N','14N','14N','14N','E','63Cu'}`.
- Lines 17: computes `inter.zeeman.eigs` using `inter.zeeman.eigs=cell(6,1)`.
- Lines 18: computes `inter.zeeman.euler` using `inter.zeeman.euler=cell(6,1)`.
- Lines 19: computes `inter.coupling.eigs` using `inter.coupling.eigs=cell(6,6)`.
- Lines 20: computes `inter.coupling.euler` using `inter.coupling.euler=cell(6,6)`.
- Lines 21: computes `inter.coupling.scalar` using `inter.coupling.scalar=cell(6,6)`.
- Lines 24: computes `inter.zeeman.eigs{5,1}` using `inter.zeeman.eigs{5,1}=[2.0509 2.0509 2.1801]`.
- Lines 25: computes `inter.zeeman.euler{5,1}` using `inter.zeeman.euler{5,1}=[0 0 0]`.
- Lines 28: computes `inter.coupling.eigs{5,6}` using `inter.coupling.eigs{5,6}=[-70.9257 -70.9257 -575.0219]*1e6`.
- Lines 29: computes `inter.coupling.euler{5,6}` using `inter.coupling.euler{5,6}=[0 0 0]`.
- Lines 30: computes `inter.coupling.scalar{1,5}` using `inter.coupling.scalar{1,5}=46.0345*1e6`.
- Lines 31: computes `inter.coupling.scalar{2,5}` using `inter.coupling.scalar{2,5}=46.0345*1e6`.
- Lines 32: computes `inter.coupling.scalar{3,5}` using `inter.coupling.scalar{3,5}=46.0345*1e6`.
- Lines 33: computes `inter.coupling.scalar{4,5}` using `inter.coupling.scalar{4,5}=46.0345*1e6`.
- Lines 36: computes `bas.formalism` using `bas.formalism='zeeman-hilb'`.
- Lines 37: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 40: computes `bas.sym_group` using `bas.sym_group={'S4'}`.

## Implementation structure

- Field swept EPR spectrum of copper porphyrin complex, computed
- by finding resonance fields and transition moments.
- Calculation time: minutes.
- Magnet field
- Isotopes
- Array preallocation
- Zeeman interactions
- Hyperfine interactions
- Basis set
- Symmetry
- Spinach housekeeping
- Experiment parameters

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `fieldsweep()`, `kfigure()`, `kxlabel()`, `kylabel()`.
