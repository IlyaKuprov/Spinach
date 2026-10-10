# examples/fundamentals/convention_tests/euler_sup_test.m

- MATLAB implementation: [examples/fundamentals/convention_tests/euler_sup_test.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/fundamentals/convention_tests/euler_sup_test.m)

- Signature: `euler_sup_test()`

## Purpose

Checks the angle-order convention and singular-branch handling of `euler_sup` by comparing its output rotation with direct direction-cosine-matrix (DCM) composition. This is a randomised consistency test, not a fit or a demonstration of every possible Euler-angle input.

## Setup and checks

Run `euler_sup_test()`. It first composes two zero triples and requires the resulting DCM to differ from the identity by no more than `1e-12` in the 1-norm. It then runs 2,000 random pairs: each Euler component is sampled as `8*pi*(rand-0.5)`, and the reference is explicitly `euler2dcm(ang_two)*euler2dcm(ang_one)`. The DCM from `euler_sup(ang_one,ang_two)` must agree with that ordered product to a 1-norm residual no greater than `1e-3`.

Three further 500-pair loops exercise near-singular branches: both middle angles are sampled near zero; both are sampled near `pi`; or the angles share their outer components and their middle angles sum to `pi`. Each uses the same ordered DCM reference and `1e-3` 1-norm threshold. These cases make the argument order and the branch-sensitive middle-angle situations particularly useful when locating a convention mismatch.

## Observable result and scope

If all checks pass, the function prints `Euler angle superposition test PASSED.`; failures raise an error at the first failed identity, random, or singular-branch comparison. It produces no plot and returns no fit or rotation result. The samples are drawn from MATLAB's current random stream; the source does not set a seed, so this is repeated randomised coverage rather than an exhaustive guarantee.
