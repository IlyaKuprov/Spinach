# examples/fundamentals/convention_tests/euler_sup_test.m

- Signature: `euler_sup_test()`

## Purpose

Tests whether `euler_sup` composes Euler-angle rotations consistently with direct direction-cosine-matrix (DCM) multiplication, including the near-singular branches.

## Checks

- Composing two zero-angle triples must give the identity DCM, to a 1-norm tolerance of 10⁻¹².
- For 2,000 random angle pairs drawn over [−4π, 4π] per angle, the DCM from `euler_sup` is compared with `euler2dcm(ang_two)*euler2dcm(ang_one)`; the 1-norm residual must not exceed 10⁻³.
- Three 500-case stress loops test rotations with both middle angles near 0, both near π, and same-phase rotations whose middle angles sum to π. Each uses the same direct-multiplication reference and 10⁻³ tolerance.
