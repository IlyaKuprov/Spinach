# examples/fundamentals/convention_tests/rotations_2.m

- MATLAB implementation: [examples/fundamentals/convention_tests/rotations_2.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/fundamentals/convention_tests/rotations_2.m)

- Signature: `rotations_2()`

## Question tested

The test compares two ways to represent the same orientation in a two-spin Hamiltonian: leave the interaction tensors and coordinates in the input frame and apply Spinach's `orientation` operator to the Hamiltonian, or rotate the tensors and coordinates in the interaction input and evaluate at zero orientation.

## Construction

It draws two shift matrices with `randn(3,3)` and a coupling matrix with `100*randn(3,3)`. The two spins are `1H` and `15N`; the magnet-field setting is `14.1`, the basis is `sphten-liouv` with approximation `none`, and the coordinates are `[0.7 0.8 0.9]` and `[1.5 2.5 3.5]`.

For construction A, it forms the lab-frame Hamiltonian from these unrotated inputs and sets `H_A = H + orientation(Q,[1 2 3])`. For construction B, it sets `R=euler2dcm(1,2,3)`, replaces each interaction matrix `A` by `R*A*R'`, rotates each coordinate row `r` as `r*R'`, then forms the Hamiltonian with `orientation(Q,[0 0 0])`. The comparison is `norm(H_A-H_B,1)`; the code errors if the residual exceeds `1e-3` and otherwise reports the residual.

The random matrices are not explicitly symmetrised in the source, so this is a convention-level comparison using generated matrices, not a validation of molecular interaction parameters. The source gives the numeric magnet-field setting but no unit annotation.
