# examples/fundamentals/derivative_tests/dirdiff_1.m

- MATLAB implementation: [examples/fundamentals/derivative_tests/dirdiff_1.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/fundamentals/derivative_tests/dirdiff_1.m)

- Signature: `dirdiff_1()`

## Purpose

Compares first and mixed second derivatives of a matrix propagator with finite differences. The source loops over the Spinach formalisms `sphten-liouv`, `zeeman-liouv`, and `zeeman-hilb`, calling `dirdiff_test_system` for each. Its comment specifies that the derivative comparison itself uses random matrices and is formalism-blind: the formalism loop exercises test-system construction rather than three distinct physical Hamiltonian cases.

## Matrix inputs and first derivative

For each constructed test system, the script draws a random complex `5×5` matrix `H` and symmetrises it as `(H+H')/2`. Two random complex direction matrices `A` and `B` are likewise symmetrised and scaled by `1/20`.

With step `1e-3` and propagator argument `1`, the first numerical derivative is the centred difference

`(propagator(spin_system,H+1e-3*A,1)-propagator(spin_system,H-1e-3*A,1))/(2e-3)`.

It is compared in the spectral norm with the second cell returned by `dirdiff(spin_system,H,A,1,2)`. The relative error is `norm(D_num-D_anl{2},2)/norm(D_num,2)`; a value above `1e-5` raises `first derivative test failed`.

## Mixed second derivative

The numerical mixed derivative is the four-corner central difference of `propagator(spin_system,H±1e-3*A±1e-3*B,1)`, divided by `4e-6`. The analytical calls are `P=dirdiff(spin_system,H,{A,B},1,3)` and `Q=dirdiff(spin_system,H,{B,A},1,3)`; the compared matrix is `(P{3}+Q{3})/2`. The relative spectral-norm error `norm(D_num-D_anl,2)/norm(D_num,2)` must not exceed `1e-3`; otherwise the script raises `second derivative test failed`.

## Scope

The matrices are randomised inside each formalism iteration and no random seed is set. The test covers one first-direction and one mixed-direction comparison per constructed system; it does not compare derivatives for a specified Spinach physical Hamiltonian.
