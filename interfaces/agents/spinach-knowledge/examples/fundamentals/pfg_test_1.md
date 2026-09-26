# examples/fundamentals/pfg_test_1.m

- Signature: `pfg_test_1()`

## Purpose

A test of the explicit gradient pulse function that uses the auxiliary matrix formalism to compute sample volume integral. For details, see:

## Physical / mathematical content

- Fundamentals examples. These are unit tests, convention checks, and pedagogical demonstrations of operator algebra, perturbation theory, tensor conventions, symmetry, quadrature, and numerical differentiation.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Implementation structure

- A test of the explicit gradient pulse function that uses the auxiliary
- matrix formalism to compute sample volume integral. For details, see:
- Magnet and isotopes
- Random chemical shifts and couplings
- Basis set
- Spinach housekeeping
- Spin Hamiltonian
- Build initial state vector
- Determine projection quantum numbers of the basis
- Determine the coherence order of each state
- Find out which coherence orders are present
- Weight coherence orders by the number of states
