# examples/quantum_tech/circuit_qed/cavity_fock_grape_a.m

Source: https://github.com/IlyaKuprov/Spinach/blob/main/examples/quantum_tech/circuit_qed/cavity_fock_grape_a.m

## Objective and physical model

This GRAPE example prepares cavity Fock state `|2>` from vacuum using a dispersively coupled qubit. A linear drive on a harmonic cavity alone produces displaced-vacuum states, not an isolated number state; the qubit's state-dependent dispersive shift supplies the nonlinearity used by the joint cavity–qubit controls. The model and parameters are attributed in the source to the bosonic GRAPE example in the para-qeet package.

The cavity is simulated in two separate truncations, with three and four Fock levels (`C3` and `C4`), each coupled to a two-level qubit (`E`). In the selected rotating frame the cavity is on resonance with its drive; the effective Hamiltonian includes the dispersive number-dependent cavity–qubit shift, set to `656.2e3` Hz (656.2 kHz). The two cavity controls are its in-phase and quadrature displacement operators, accompanied by qubit `Lx` and `Ly` controls. The initial state is cavity vacuum with the qubit in its upper state; the target is cavity `|2>` with the qubit still upper. The calculation is a closed-system coherent model: the source specifies no relaxation, dephasing, or measurement process.

## Controls and truncation comparison

Four piecewise-constant controls act on the cavity's two quadratures and the qubit's x and y operators. Each pulse has 40 slices of 33 ns, for a total duration of 1.32 μs. The optimisation uses L-BFGS, at most 300 iterations, the source's `NS` penalty with weight 0.001, and control scaling `1.76828e7` (the source does not annotate a physical unit for this scale). The initial guess has a Gaussian envelope centred at 0.66 μs, with 1.32 μs divided by 8 as its width parameter; its in-phase cavity and qubit channels start at amplitude 0.7 and the quadrature channels at zero.

The script optimises separately in `C3` and `C4`, directly propagates each resulting pulse to calculate transfer fidelity, and evaluates the `C3`-optimised pulse again in `C4`. It requires each in-space optimisation fidelity to reach 0.95, then reports the difference between the `C4`-optimised fidelity and the transplanted `C3` pulse's fidelity. This cross-test probes sensitivity to the chosen Fock truncation; two finite truncations alone do not prove convergence. The source comments describe the expected underperformance of the transplanted pulse, but no MATLAB result is included here and no numerical fidelity or convergence claim is made.
