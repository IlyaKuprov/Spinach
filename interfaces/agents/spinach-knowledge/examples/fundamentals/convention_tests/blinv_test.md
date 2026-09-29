# examples/fundamentals/convention_tests/blinv_test.m

- MATLAB implementation: [examples/fundamentals/convention_tests/blinv_test.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/fundamentals/convention_tests/blinv_test.m)

- Signature: `blinv_test()`

## Purpose

Checks second-rank Blicharski self- and cross-invariants against products of the five spherical-tensor coefficients returned by `mat2sphten`. It is an internal convention/consistency test on randomly drawn 3×3 matrices.

## Convention checked

The function draws `A=rand(3,3)` and `B=rand(3,3)`, takes the third output of `mat2sphten` as `Phi_A` and `Phi_B`, obtains each self-invariant `Dsq` from `blinv`, and obtains the cross-invariant `X_AB` from `blprod`. For each tensor, the rank-2 coefficient combination

`Phi(3)^2-2*Phi(2)*Phi(4)+2*Phi(1)*Phi(5)`

is compared with `(2/3)*Dsq`. For the pair, the combination

`Phi_A(1)*Phi_B(5)-Phi_A(2)*Phi_B(4)+Phi_A(3)*Phi_B(3)-Phi_A(4)*Phi_B(2)+Phi_A(5)*Phi_B(1)`

is compared with `(2/3)*X_AB`. These explicit index/sign patterns are the convention exercised by the test.

## Pass criterion and scope

The two self residuals are tested together, then the cross residual is tested separately; each test fails if its absolute residual is greater than `10*eps`, otherwise it prints the corresponding passed message. There are no inputs or returned values, and the random-number generator is not seeded. Since the matrices are sampled by unconstrained `rand(3,3)`, this example checks these algebraic identities on random matrices; it does not separately impose or test symmetry, tracelessness, or any physical tensor constraints.
