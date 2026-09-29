# examples/fundamentals/state_tests/state_consistency_2.m

Source: [examples/fundamentals/state_tests/state_consistency_2.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/fundamentals/state_tests/state_consistency_2.m)

Signature: `state_consistency_2()`

## Tested question

Does one constructed two-spin operator give the same Hilbert-space matrix when built directly in Zeeman Hilbert space, built in Zeeman Liouville space, or built in spherical-tensor Liouville space and projected back to the Zeeman basis?

## System and state

The source sets `sys.magnet=14.1`, isotopes `14N` and `235U`, and both scalar Zeeman entries to zero. It declares no couplings. The basis approximation is `none`. The test state is `state(...,{'Lz','Lx'},{1,2}) + state(...,{'L+'},{1})`: the product component has `Lz` on spin 1 and `Lx` on spin 2, with an additional `L+` component on spin 1.

## Representation comparison

The source builds that state in `zeeman-hilb`, `zeeman-liouv`, and `sphten-liouv`. The Zeeman-Liouville vector is reshaped to a 24-by-24 matrix. The spherical-tensor vector is first multiplied by `sphten2zeeman(spin_system)`, then reshaped to the same dimensions. The three Hilbert-space matrices are compared pairwise using the matrix 1-norm; any difference greater than `1e-6` triggers the failure branch.

## Output and scope

The function prints either `State consistency test PASSED.` or raises `State consistency test FAILED.`; it does not return the matrices or report their residuals. This is a fixed state-space representation check, not a quadrature test, and does not establish agreement for other systems, states, or approximations.
