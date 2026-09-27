# examples/quantum_tech/circuit_qed/cavity_fock_grape_a.m

- Signature: `cavity_fock_grape_a()`

## Purpose

Uses GRAPE to prepare cavity Fock state 2 from the cavity vacuum through a dispersively coupled qubit. It optimises and compares pulses at cavity truncations `C3` and `C4`, then tests the smaller-space pulse in the larger space. The source attributes the model and parameters to the bosonic GRAPE example in the paraqeet package and gives a calculation time of minutes.

## Physical / mathematical content

A linear drive alone cannot make a Fock state from the vacuum of a harmonic cavity; the dispersive coupling to the qubit supplies the required nonlinearity. The model uses a 656.2 kHz dispersive coupling. The initial state is cavity vacuum with the qubit in its upper level; the target is cavity Fock state 2 with the qubit in the same level. The controls are the two cavity quadratures and the qubit `Lx` and `Ly` operators. Comparing the two truncations demonstrates the source's point that a pulse optimised in a smaller Fock space can lose transfer fidelity when evaluated in a larger one.

## Numerical / algorithmic content

For each truncation, the script optimises a 40-slice pulse with 33 ns per slice using `fmaxnewton` and `grape_xy`. The control settings are power level `1.76828e7`, the `NS` penalty with weight 0.001, L-BFGS, and a maximum of 300 iterations. A Gaussian initial guess is applied to the in-phase cavity and qubit channels. Transfer fidelity is recomputed by direct slice-by-slice propagation; both optimised pulses must reach at least 0.95 fidelity.

## Implementation structure

The script constructs each Zeeman Hilbert-space model and its cavity-frame drift Hamiltonian, defines the four control operators and initial/target states, runs the two optimisations, and evaluates all three comparisons: optimise and test at `C3`, optimise and test at `C4`, and optimise at `C3` but test at `C4`.
