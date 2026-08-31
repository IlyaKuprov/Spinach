# examples/fundamentals/derivative_tests/dirdiff_1.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/fundamentals/derivative_tests/dirdiff_1.m`
- Signature: `dirdiff_1()`
- Total lines: 58

## Purpose

Test of matrix exponential differentiation routines. Analytical derivatives are compared to central finite differences.

## Physical / mathematical content

- Fundamentals examples. These are unit tests, convention checks, and pedagogical demonstrations of operator algebra, perturbation theory, tensor conventions, symmetry, quadrature, and numerical differentiation.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 8-9: Formalisms to test; implemented by `formalisms={'sphten-liouv','zeeman-liouv','zeeman-hilb'}`.
- Lines 11-12: Loop over formalisms; implemented by `for n=1:numel(formalisms)`.
- Lines 14-15: Get the Spinach object; implemented by `spin_system=dirdiff_test_system(formalisms{n})`.
- Lines 17-18: Random Hamiltonian; implemented by `H=randn(5)+1i*randn(5); H=(H+H')/2`.
- Lines 20-21: Random direction operators; implemented by `A=randn(5)+1i*randn(5); A=(A+A')/20`.
- Lines 24-26: First derivative, numerical; implemented by `D_num=(propagator(spin_system,H+1e-3*A,1)- propagator(spin_system,H-1e-3*A,1))/2e-3`.
- Lines 28-29: First derivative, analytical; implemented by `D_anl=dirdiff(spin_system,H,A,1,2)`.
- Lines 31-32: Test the first derivative; implemented by `if norm(D_num-D_anl{2},2)/norm(D_num,2)>1e-5`.
- Lines 38-42: Second derivative, numerical; implemented by `D_num=(propagator(spin_system,H+1e-3*A+1e-3*B,1)- propagator(spin_system,H+1e-3*A-1e-3*B,1)- propagator(spin_system,H-1e-3*A+1e-3*B,1)+ propagator(spin_system,H-1e-3*A-1…`.
- Lines 44-45: Second derivative, analytical; implemented by `P=dirdiff(spin_system,H,{A,B},1,3)`.
- Lines 49-50: Test the second derivative; implemented by `if norm(D_num-D_anl,2)/norm(D_num,2)>1e-3`.

### Control flow inferred from the code

- Line 12: `for` loop over `n=1:numel(formalisms)`.
- Line 32: conditional branch on `norm(D_num-D_anl{2},2)/norm(D_num,2)>1e-5`.
- Line 50: conditional branch on `norm(D_num-D_anl,2)/norm(D_num,2)>1e-3`.

### Key state/data transformations

- Lines 9: computes `formalisms` using `formalisms={'sphten-liouv','zeeman-liouv','zeeman-hilb'}`.
- Lines 15: computes `spin_system` using `spin_system=dirdiff_test_system(formalisms{n})`.
- Lines 18: computes `H` using `H=randn(5)+1i*randn(5); H=(H+H')/2`.
- Lines 21: computes `A` using `A=randn(5)+1i*randn(5); A=(A+A')/20`.
- Lines 22: computes `B` using `B=randn(5)+1i*randn(5); B=(B+B')/20`.
- Lines 25-26: computes `D_num` using `D_num=(propagator(spin_system,H+1e-3*A,1)- propagator(spin_system,H-1e-3*A,1))/2e-3`.
- Lines 29: computes `D_anl` using `D_anl=dirdiff(spin_system,H,A,1,2)`.
- Lines 45: computes `P` using `P=dirdiff(spin_system,H,{A,B},1,3)`.
- Lines 46: computes `Q` using `Q=dirdiff(spin_system,H,{B,A},1,3)`.

## Implementation structure

- Test of matrix exponential differentiation routines. Analytical
- derivatives are compared to central finite differences.
- Formalisms to test
- Loop over formalisms
- Get the Spinach object
- Random Hamiltonian
- Random direction operators
- First derivative, numerical
- First derivative, analytical
- Test the first derivative
- Second derivative, numerical
- Second derivative, analytical

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `dirdiff_test_system()`, `propagator()`, `dirdiff()`.
