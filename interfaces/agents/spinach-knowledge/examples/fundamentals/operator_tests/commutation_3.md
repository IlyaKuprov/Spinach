# examples/fundamentals/operator_tests/commutation_3.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/fundamentals/operator_tests/commutation_3.m`
- Signature: `commutation_3()`
- Total lines: 51

## Purpose

Commutators of simple operators and superoperators. The test calculation is performed three times in the three formalisms supported by Spinach.

## Physical / mathematical content

- Fundamentals examples. These are unit tests, convention checks, and pedagogical demonstrations of operator algebra, perturbation theory, tensor conventions, symmetry, quadrature, and numerical differentiation.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 9-10: System specification; implemented by `sys.magnet=14.1`.
- Lines 16-17: Preallocate the answer; implemented by `answer=zeros(7,3,'like',1i)`.
- Lines 19-20: Run the tests; implemented by `formalisms={'zeeman-hilb','zeeman-liouv','sphten-liouv'}`.
- Lines 43-44: Report the outcome; implemented by `if norm(answer,'fro')<1e-6`.

### Control flow inferred from the code

- Line 21: `for` loop over `n=1:numel(formalisms)`.
- Line 44: conditional branch on `norm(answer,'fro')<1e-6`.

### Key state/data transformations

- Lines 10: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 11: computes `sys.isotopes` using `sys.isotopes={'1H','235U'}`.
- Lines 12: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={2.5 1.0}`.
- Lines 13: computes `inter.coupling.scalar{1,2}` using `inter.coupling.scalar{1,2}=10`.
- Lines 14: computes `inter.coupling.scalar{2,1}` using `inter.coupling.scalar{2,1}=10`.
- Lines 17: computes `answer` using `answer=zeros(7,3,'like',1i)`.
- Lines 20: computes `formalisms` using `formalisms={'zeeman-hilb','zeeman-liouv','sphten-liouv'}`.
- Lines 22: computes `bas.formalism` using `bas.formalism=formalisms{n}`.
- Lines 23: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 24: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 26: computes `Up` using `Up=operator(spin_system,'L+','235U')`.
- Lines 27: computes `Um` using `Um=operator(spin_system,'L-','235U')`.
- Lines 28: computes `Ux` using `Ux=operator(spin_system,'Lx','235U')`.
- Lines 29: computes `Uy` using `Uy=operator(spin_system,'Ly','235U')`.
- Lines 30: computes `Uz` using `Uz=operator(spin_system,'Lz','235U')`.
- Lines 31: computes `HpUp` using `HpUp=operator(spin_system,{'L+','L+'},{1,2})`.
- Lines 32: computes `HmUm` using `HmUm=operator(spin_system,{'L-','L-'},{1,2})`.
- Lines 33: computes `Hz` using `Hz=operator(spin_system,'Lz','1H')`.

## Implementation structure

- Commutators of simple operators and superoperators. The test
- calculation is performed three times in the three formalisms
- supported by Spinach.
- System specification
- Preallocate the answer
- Run the tests
- Report the outcome

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `operator()`, `answer()`.
