# examples/quantum_tech/circuit_qed/cavity_fock_grape_b.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/quantum_tech/circuit_qed/cavity_fock_grape_b.m`
- Signature: `cavity_fock_grape_b()`
- Total lines: 113

## Purpose

GRAPE preparation of a cavity Fock state through a dispersively coupled qubit using smooth band-limited drives. The controls are expanded in an orthonormal basis of slow sine and cosine waves, and the optimisation runs over the expansion coefficients, so that the resulting pulses are hardware-friendly smooth envelopes rather than free piecewise-constant switches. Compare with the piecewise-constant treatment in cavit

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
- coupled qubit using smooth band-limited drives. The controls are
- expanded in an orthonormal basis of slow sine and cosine waves,
- and the optimisation runs over the expansion coefficients, so
- that the resulting pulses are hardware-friendly smooth envelopes
- rather than free piecewise-constant switches. Compare with the
- piecewise-constant treatment in cavity_fock_grape_a.m; model and
- parameters follow the smooth pulse bosonic GRAPE example of the
- paraqeet package.
- Calculation time: minutes
- Magnet field
- Truncated cavity mode and a qubit

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `hamiltonian()`, `assume()`, `operator()`, `state()`, `wave_basis()`, `optimcon()`, `fmaxnewton()`, `pulse()`, `propagator()`, `cumsum()`, `kfigure()`, `kxlabel()`, `kylabel()`, `ktitle()`.
