# examples/fundamentals/state_tests/state_consistency_3.m

Source: [examples/fundamentals/state_tests/state_consistency_3.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/fundamentals/state_tests/state_consistency_3.m)

Signature: `state_consistency_3()`

## Tested question

Do Spinach's deuterium-pair singlet, triplet, and quintet projectors match explicit spin-1 component-state projectors, and do all nine projectors sum to the two-spin unit state?

## System and reference states

The source uses two `2H` spins, zero magnet field, and zero scalar Zeeman entries, with basis approximation `none` and `zeeman-hilb` formalism. `deut_pair(spin_system,1,2)` supplies one singlet projector, three triplet projectors, and five quintet projectors. For the independent references, the source defines the three single-spin components `alp`, `bet`, and `gam`, forms the singlet, triplet, and quintet component vectors, and makes each projector as a vector times its conjugate transpose. The component definitions are cited to Eq. 1 of [the 1998 paper](https://doi.org/10.1016/S0009-2614(98)00784-2).

## Checks and output

Each Spinach projector is compared with its corresponding explicit projector using the matrix 1-norm and threshold `1e-6`. The source then rebuilds the same spin system in Zeeman Hilbert space, obtains the projectors again, and compares their sum with `state(spin_system,{'E','E'},{1 2})`, the two-spin unit state, using the same norm and threshold. A larger residual raises `State construction test FAILED.`; otherwise the function prints `State construction test PASSED.`

## Scope

This checks projector construction and completeness for this zero-field spin-1 pair; it does not test quadrature or establish state identities for other spin pairs or Hamiltonians.
