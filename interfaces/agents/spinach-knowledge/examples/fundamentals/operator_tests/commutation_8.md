# examples/fundamentals/operator_tests/commutation_8.m

- Signature: `commutation_8()`

## Purpose

Tests Hilbert-to-Liouville operator actions and the first-rank Stevens-to-Pauli mapping.

## Physical / mathematical content

For random complex 6-by-6 matrices `H` and `R`, the script checks that the left, right, commutator, and anticommutator Liouville operators generated from `H` act on the vectorized `R` as `H*R`, `R*H`, `H*R-R*H`, and `H*R+R*H`. It also compares first-rank Stevens operators with the Pauli Cartesian operators at multiplicities 2, 3, 5, and 8.

## Numerical / algorithmic content

All action and operator-mapping residuals use the fixed threshold `1e-10`. Liouville action residuals are measured with the vector 2-norm; Stevens/Pauli differences use the Frobenius norm. A residual above threshold raises an error.

## Implementation structure

The first test builds the four Liouville representations and compares their action on `hilb2liouv(R,'statevec')` with explicitly formed Hilbert-space products. The second loops over multiplicities, obtains the Pauli and rank-1 Stevens operators, and checks `O_10=L.z`, `O_11=L.x`, and `O_1m1=L.y` before reporting success.
