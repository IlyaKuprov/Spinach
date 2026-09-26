# examples/quantum_tech/circuit_qed/cavity_fock_grape_a.m

- Signature: `cavity_fock_grape_a()`

## Purpose

GRAPE preparation of a cavity Fock state through a dispersively coupled qubit, using piecewise-constant drives on both the cavity and the qubit. A linear drive alone cannot make a Fock state out of the vacuum of a harmonic mode; the qubit conditions the cavity phase through the dispersive shift and thereby provides the requi- red nonlinearity. The optimisation is run at two Fock space trun- cations; the pulse optimis

## Physical / mathematical content

- Quantum-technology examples. The files in this area model cavity QED, transmon qubits, NV centres, and related effective Hamiltonians. The recurring mathematics is finite-dimensional quantum dynamics with ladder operators, rotating-wave-style couplings, anharmonic oscillator terms, avoided crossings, and coherent control in coupled few-mode systems.
- The numerical method is quasi-Newton optimisation: curvature information is approximated from successive step and gradient differences instead of forming exact second derivatives every iteration.
- The numerical method is limited-memory quasi-Newton optimisation, which keeps only a short curvature history and is therefore suitable for waveform vectors too large for dense Hessians.
- The optimisation logic is Newton or Newton-like: search directions use first- and second-order local curvature information, usually with regularisation or line-search safeguards.
- The control theory content is GRAPE: fidelity derivatives are propagated through a piecewise-constant pulse sequence so that waveform samples can be improved by gradient-based optimisation.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.

## Implementation structure

- GRAPE preparation of a cavity Fock state through a dispersively
- coupled qubit, using piecewise-constant drives on both the cavity
- and the qubit. A linear drive alone cannot make a Fock state out
- of the vacuum of a harmonic mode; the qubit conditions the cavity
- phase through the dispersive shift and thereby provides the requi-
- red nonlinearity. The optimisation is run at two Fock space trun-
- cations; the pulse optimised in the smaller space underperforms
- when it is re-evaluated in the larger one -optimal control solu-
- tions must be converged with respect to the Fock space truncation.
- Model and parameters from the bosonic GRAPE example of the para-
- qeet package.
- Calculation time: minutes
