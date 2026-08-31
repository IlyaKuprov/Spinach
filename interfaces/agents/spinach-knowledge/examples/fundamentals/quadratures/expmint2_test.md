# examples/fundamentals/quadratures/expmint2_test.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/fundamentals/quadratures/expmint2_test.m`
- Signature: `expmint2_test()`
- Total lines: 48

## Purpose

Verification of expmint2 against numerical integration.

## Physical / mathematical content

- Fundamentals examples. These are unit tests, convention checks, and pedagogical demonstrations of operator algebra, perturbation theory, tensor conventions, symmetry, quadrature, and numerical differentiation.
- Orientation or trajectory averaging is performed numerically, so grid design, weights, and integration error control matter directly to accuracy and runtime.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 9-10: Random dimension and upper limit; implemented by `n=10+randi(5); ul=randi(10)+rand()`.
- Lines 12-13: Generate random matrices; implemented by `A=rand(n)+1i*rand(n); A=(A+A')/2`.
- Lines 19-20: Bootstrap Spinach; implemented by `spin_system=bootstrap()`.
- Lines 22-23: Call Spinach function; implemented by `int_spinach=expmint2(spin_system,A,B,C,D,E,ul)`.
- Lines 25-26: Inner integrand and its integral; implemented by `fun_inner=@(x)expm(1i*C*x)*D*expm(-1i*E*x)`.
- Lines 29-30: Outer integrand and its integral; implemented by `fun_outer=@(t)expm(1i*A*t)*B*int_inner(t)`.
- Lines 33-34: Call Matlab integrator; implemented by `int_matlab=int_outer(ul)`.
- Lines 36-37: Test the error norm; implemented by `norm_diff=norm(int_spinach-int_matlab,'fro')`.

### Control flow inferred from the code

- Line 39: conditional branch on `norm_diff/norm_comp>10*n*eps('double')`.

### Key state/data transformations

- Lines 10: computes `n` using `n=10+randi(5); ul=randi(10)+rand()`.
- Lines 13: computes `A` using `A=rand(n)+1i*rand(n); A=(A+A')/2`.
- Lines 14: computes `C` using `C=rand(n)+1i*rand(n); C=(C+C')/2`.
- Lines 15: computes `E` using `E=rand(n)+1i*rand(n); E=(E+E')/2`.
- Lines 16: computes `B` using `B=rand(n)+1i*rand(n)`.
- Lines 17: computes `D` using `D=rand(n)+1i*rand(n)`.
- Lines 20: computes `spin_system` using `spin_system=bootstrap()`.
- Lines 23: computes `int_spinach` using `int_spinach=expmint2(spin_system,A,B,C,D,E,ul)`.
- Lines 26: computes `fun_inner` using `fun_inner=@(x)expm(1i*C*x)*D*expm(-1i*E*x)`.
- Lines 27: computes `int_inner` using `int_inner=@(t)expm(-1i*C*t)*integral(fun_inner,0,t,'ArrayValued',true)`.
- Lines 30: computes `fun_outer` using `fun_outer=@(t)expm(1i*A*t)*B*int_inner(t)`.
- Lines 31: computes `int_outer` using `int_outer=@(T)expm(-1i*A*T)*integral(fun_outer,0,T,'ArrayValued',true)`.
- Lines 34: computes `int_matlab` using `int_matlab=int_outer(ul)`.
- Lines 37: computes `norm_diff` using `norm_diff=norm(int_spinach-int_matlab,'fro')`.
- Lines 38: computes `norm_comp` using `norm_comp=max([norm(int_spinach,'fro') norm(int_spinach,'fro')])`.

## Implementation structure

- Verification of expmint2 against numerical
- integration.
- Random dimension and upper limit
- Generate random matrices
- Bootstrap Spinach
- Call Spinach function
- Inner integrand and its integral
- Outer integrand and its integral
- Call Matlab integrator
- Test the error norm

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `randi()`, `bootstrap()`, `expmint2()`, `integral()`, `int_inner()`, `int_outer()`, `eps()`.
