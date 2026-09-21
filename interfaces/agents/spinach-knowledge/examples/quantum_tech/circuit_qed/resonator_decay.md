# examples/quantum_tech/circuit_qed/resonator_decay.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/quantum_tech/circuit_qed/resonator_decay.m`
- Signature: `resonator_decay()`
- Total lines: 117

## Purpose

Open-system dynamics of a leaky microwave resonator at finite temperature. A Fock state decays as a downward cascade through the level ladder, and a coherent state decays with its Poisson population structure largely preserved; both settle into the thermal state of the mode. The mean photon number follows the analytical amplitude damping solution in both cases, up to the distortion of the weak thermal channel by the 

## Physical / mathematical content

- Quantum-technology examples. The files in this area model cavity QED, transmon qubits, NV centres, and related effective Hamiltonians. The recurring mathematics is finite-dimensional quantum dynamics with ladder operators, rotating-wave-style couplings, anharmonic oscillator terms, avoided crossings, and coherent control in coupled few-mode systems.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

## Implementation structure

- Open-system dynamics of a leaky microwave resonator at finite
- temperature. A Fock state decays as a downward cascade through
- the level ladder, and a coherent state decays with its Poisson
- population structure largely preserved; both settle into the
- thermal state of the mode. The mean photon number follows the
- analytical amplitude damping solution in both cases, up to the
- distortion of the weak thermal channel by the Fock space trun-
- cation. Model and parameters from the resonator decay example
- of the paraqeet package.
- Calculation time: seconds
- Magnet field
- Microwave resonator with five Fock levels

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `hamiltonian()`, `assume()`, `relaxation()`, `state()`, `coherent()`, `pops_fock()`, `pops_coh()`, `n_fock()`, `n_coh()`, `kfigure()`, `scale_figure()`, `subplot()`, `kxlabel()`, `kylabel()`.
