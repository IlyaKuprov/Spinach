# tests/kernel/test_scalar_coupling_hamiltonian.m

- Signature: `result=test_scalar_coupling_hamiltonian()`

## Purpose

Checks the isotropic scalar-coupling Hamiltonian for two protons with a 10 Hz coupling.

## Physical / mathematical content

For isotropic coupling J, the test checks the rotationally invariant Hamiltonian `2*pi*J*(Ix*Sx+Iy*Sy+Iz*Sz)`, expressed in rad/s.

## Numerical / algorithmic content

Builds the two-proton Zeeman-Hilbert spin system with zero chemical shifts and a 10 Hz scalar coupling. It constructs the Hamiltonian under the NMR assumption, forms the explicit spin-operator reference, and compares the matrices using absolute and relative tolerances of `1e-9` and `1e-12`.

## Outputs

`result` is the regression-test record, including the comparison result and explanatory messages.

## Implementation structure

The test builds the spin system and coupling operators, forms the textbook reference Hamiltonian, and checks it against Spinach's Hamiltonian.
