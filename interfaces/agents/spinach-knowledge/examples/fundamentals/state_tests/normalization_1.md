# examples/fundamentals/state_tests/normalization_1.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/fundamentals/state_tests/normalization_1.m`
- Signature: `normalization_1()`
- Total lines: 47

## Purpose

Internal consistency test for the state vectors and matrices in each of the three formalisms supported by Spinach.

## Physical / mathematical content

- Fundamentals examples. These are unit tests, convention checks, and pedagogical demonstrations of operator algebra, perturbation theory, tensor conventions, symmetry, quadrature, and numerical differentiation.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 8-9: System specification; implemented by `sys.magnet=14.1`.
- Lines 16-17: Preallocate the answer; implemented by `norm_diffs=zeros(6,3,'like',1i)`.
- Lines 19-20: Compute norm differences; implemented by `formalisms={'zeeman-hilb','zeeman-liouv','sphten-liouv'}`.
- Lines 39-40: Display the answers; implemented by `if any(abs(norm_diffs)>1e-6,'all')`.

### Control flow inferred from the code

- Line 21: `for` loop over `n=1:numel(formalisms)`.
- Line 40: conditional branch on `any(abs(norm_diffs)>1e-6,'all')`.

### Key state/data transformations

- Lines 9: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 10: computes `sys.isotopes` using `sys.isotopes={'1H','235U'}`.
- Lines 11: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={2.5 1.0}`.
- Lines 12: computes `inter.coupling.scalar{1,2}` using `inter.coupling.scalar{1,2}=10`.
- Lines 13: computes `inter.coupling.scalar{2,1}` using `inter.coupling.scalar{2,1}=10`.
- Lines 14: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 17: computes `norm_diffs` using `norm_diffs=zeros(6,3,'like',1i)`.
- Lines 20: computes `formalisms` using `formalisms={'zeeman-hilb','zeeman-liouv','sphten-liouv'}`.
- Lines 22: computes `bas.formalism` using `bas.formalism=formalisms{n}`.
- Lines 23: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 25: computes `Ux` using `Ux=state(spin_system,'Lx','235U')`.
- Lines 26: computes `Uy` using `Uy=state(spin_system,'Ly','235U')`.
- Lines 27: computes `Uz` using `Uz=state(spin_system,'Lz','235U')`.
- Lines 28: computes `Hx` using `Hx=state(spin_system,'Lx','1H')`.
- Lines 29: computes `Hy` using `Hy=state(spin_system,'Ly','1H')`.
- Lines 30: computes `Hz` using `Hz=state(spin_system,'Lz','1H')`.
- Lines 31: computes `norm_diffs(1,n)` using `norm_diffs(1,n)=norm(full(Ux))-norm(full(Uy))`.
- Lines 32: computes `norm_diffs(2,n)` using `norm_diffs(2,n)=norm(full(Uy))-norm(full(Uz))`.

## Implementation structure

- Internal consistency test for the state vectors and matrices
- in each of the three formalisms supported by Spinach.
- System specification
- Preallocate the answer
- Compute norm differences
- Display the answers

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `norm_diffs()`, `any()`.
