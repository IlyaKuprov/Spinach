# examples/fundamentals/state_tests/normalization_2.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/fundamentals/state_tests/normalization_2.m`
- Signature: `normalization_2()`
- Total lines: 52

## Purpose

Internal consistency test for the state vectors and matrices. Checks that the inner products are consistent between the for- malisms supported by Spinach. Output should be a 6x3 matrix with identical columns.

## Physical / mathematical content

- Fundamentals examples. These are unit tests, convention checks, and pedagogical demonstrations of operator algebra, perturbation theory, tensor conventions, symmetry, quadrature, and numerical differentiation.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 11-12: System specification; implemented by `sys.magnet=14.1`.
- Lines 19-20: Preallocate the answer; implemented by `norms=zeros(6,3,'like',1i)`.
- Lines 22-23: Get the norms; implemented by `formalisms={'zeeman-hilb','zeeman-liouv','sphten-liouv'}`.
- Lines 42-43: Run the tests; implemented by `mult_factor=prod(spin_system.comp.mults)`.

### Control flow inferred from the code

- Line 24: `for` loop over `n=1:numel(formalisms)`.
- Line 44: conditional branch on `norm(norms(:,1)-norms(:,2),1)>1e-6||`.

### Key state/data transformations

- Lines 12: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 13: computes `sys.isotopes` using `sys.isotopes={'1H','235U'}`.
- Lines 14: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={2.5 1.0}`.
- Lines 15: computes `inter.coupling.scalar{1,2}` using `inter.coupling.scalar{1,2}=10`.
- Lines 16: computes `inter.coupling.scalar{2,1}` using `inter.coupling.scalar{2,1}=10`.
- Lines 17: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 20: computes `norms` using `norms=zeros(6,3,'like',1i)`.
- Lines 23: computes `formalisms` using `formalisms={'zeeman-hilb','zeeman-liouv','sphten-liouv'}`.
- Lines 25: computes `bas.formalism` using `bas.formalism=formalisms{n}`.
- Lines 26: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 28: computes `Ux` using `Ux=state(spin_system,'Lx','235U')`.
- Lines 29: computes `Uy` using `Uy=state(spin_system,'Ly','235U')`.
- Lines 30: computes `Uz` using `Uz=state(spin_system,'Lz','235U')`.
- Lines 31: computes `Hx` using `Hx=state(spin_system,'Lx','1H')`.
- Lines 32: computes `Hy` using `Hy=state(spin_system,'Ly','1H')`.
- Lines 33: computes `Hz` using `Hz=state(spin_system,'Lz','1H')`.
- Lines 34: computes `norms(1,n)` using `norms(1,n)=trace(full(Ux)'*full(Ux))`.
- Lines 35: computes `norms(2,n)` using `norms(2,n)=trace(full(Uy)'*full(Uy))`.

## Implementation structure

- Internal consistency test for the state vectors and matrices.
- Checks that the inner products are consistent between the for-
- malisms supported by Spinach. Output should be a 6x3 matrix
- with identical columns.
- System specification
- Preallocate the answer
- Get the norms
- Run the tests

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `norms()`.
