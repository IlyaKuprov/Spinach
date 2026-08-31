# examples/fundamentals/operator_tests/commutation_1.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/fundamentals/operator_tests/commutation_1.m`
- Signature: `commutation_1()`
- Total lines: 42

## Purpose

Commutators of simple operators and superoperators. The test calculation is performed three times in the three formalisms supported by Spinach.

## Physical / mathematical content

- Fundamentals examples. These are unit tests, convention checks, and pedagogical demonstrations of operator algebra, perturbation theory, tensor conventions, symmetry, quadrature, and numerical differentiation.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 9-10: Simple 1-spin system; implemented by `sys.magnet=0`.
- Lines 14-15: Preallocate the answer; implemented by `answer=zeros(3,3,'like',1i)`.
- Lines 17-18: Run the tests; implemented by `formalisms={'zeeman-hilb','zeeman-liouv','sphten-liouv'}`.
- Lines 34-35: Report the outcome; implemented by `if norm(answer,'fro')<1e-6`.

### Control flow inferred from the code

- Line 19: `for` loop over `n=1:numel(formalisms)`.
- Line 35: conditional branch on `norm(answer,'fro')<1e-6`.

### Key state/data transformations

- Lines 10: computes `sys.magnet` using `sys.magnet=0`.
- Lines 11: computes `sys.isotopes` using `sys.isotopes={'1H'}`.
- Lines 12: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={0}`.
- Lines 15: computes `answer` using `answer=zeros(3,3,'like',1i)`.
- Lines 18: computes `formalisms` using `formalisms={'zeeman-hilb','zeeman-liouv','sphten-liouv'}`.
- Lines 20: computes `bas.formalism` using `bas.formalism=formalisms{n}`.
- Lines 21: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 22: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 24: computes `Lx` using `Lx=operator(spin_system,'Lx','1H')`.
- Lines 25: computes `Ly` using `Ly=operator(spin_system,'Ly','1H')`.
- Lines 26: computes `Lz` using `Lz=operator(spin_system,'Lz','1H')`.
- Lines 27: computes `Lp` using `Lp=operator(spin_system,'L+','1H')`.
- Lines 28: computes `Lm` using `Lm=operator(spin_system,'L-','1H')`.
- Lines 29: computes `answer(1,n)` using `answer(1,n)=norm(Lz*Lp-Lp*Lz-Lp,'fro')`.
- Lines 30: computes `answer(2,n)` using `answer(2,n)=norm(Lz*Lm-Lm*Lz+Lm,'fro')`.
- Lines 31: computes `answer(3,n)` using `answer(3,n)=norm(Lx*Ly-Ly*Lx-1i*Lz,'fro')`.

## Implementation structure

- Commutators of simple operators and superoperators. The test
- calculation is performed three times in the three formalisms
- supported by Spinach.
- Simple 1-spin system
- Preallocate the answer
- Run the tests
- Report the outcome

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `operator()`, `answer()`.
