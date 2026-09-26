# examples/quantum_tech/circuit_qed/cavity_dual_rail_safe.m

- Signature: `cavity_dual_rail_safe()`

## Purpose

Flux-noise dephasing rates of two dual-rail qubits whose flux-tuna- ble transmon-coupled rails share one transmon ancilla, and their suppression by a Stark-assisted flux-noise evasion (SAFE) drive on the transmon, Sec. 4.4.2 and Fig. 4.4(c) of Yunwei Lu's PhD thesis (Northwestern University, 2026). The logical dephasing rate of each dual-rail qubit is set by the flux sensitivity of the dressed sin- gle-photon transition frequency of its transmon-coupled rail, Eqs. (4.93) and (4.95); the rates are computed from the eigenvalues of the rotating frame Hamiltonian of Eq. (4.16) at a fixed drive amp- litude as functions of the transmon-drive detuning. Both rates go through a minimum in the same detuning window, so that one drive protects both qubits. Calculation time: seconds

## Physical / mathematical content

- Quantum-technology examples. The files in this area model cavity QED, transmon qubits, NV centres, and related effective Hamiltonians. The recurring mathematics is finite-dimensional quantum dynamics with ladder operators, rotating-wave-style couplings, anharmonic oscillator terms, avoided crossings, and coherent control in coupled few-mode systems.
- The effective hardware model is a weakly anharmonic oscillator. Duffing nonlinearity breaks equal level spacing and allows qubit-like addressability within a truncated bosonic ladder.

## Numerical / algorithmic content

- An eigenvalue problem is solved or analysed, so the file is extracting spectra, stationary states, avoided crossings, or modal structure from the effective Hamiltonian or superoperator.

## Implementation structure

- Flux-noise dephasing rates of two dual-rail qubits whose flux-tuna-
- ble transmon-coupled rails share one transmon ancilla, and their
- suppression by a Stark-assisted flux-noise evasion (SAFE) drive on
- the transmon, Sec. 4.4.2 and Fig. 4.4(c) of Yunwei Lu's PhD thesis
- (Northwestern University, 2026). The logical dephasing rate of each
- dual-rail qubit is set by the flux sensitivity of the dressed sin-
- gle-photon transition frequency of its transmon-coupled rail, Eqs.
- (4.93) and (4.95); the rates are computed from the eigenvalues of
- the rotating frame Hamiltonian of Eq. (4.16) at a fixed drive amp-
- litude as functions of the transmon-drive detuning. Both rates go
- through a minimum in the same detuning window, so that one drive
- protects both qubits.
