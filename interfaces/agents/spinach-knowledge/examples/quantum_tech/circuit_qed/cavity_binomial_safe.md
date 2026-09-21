# examples/quantum_tech/circuit_qed/cavity_binomial_safe.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/quantum_tech/circuit_qed/cavity_binomial_safe.m`
- Signature: `cavity_binomial_safe()`
- Total lines: 232

## Purpose

Binomial bosonic code |0L>=(|0>+|4>)/sqrt(2), |1L>=|2> in a cavity dispersively coupled to a flux-tunable transmon ancilla, and the protection of its coherences from 1/f flux noise by a Stark-assis- ted flux-noise evasion (SAFE) drive on the transmon, Sec. 4.4.1 and Fig. 4.4(a,b) of Yunwei Lu's PhD thesis (Northwestern University, 2026). The flux noise dephasing rates of the code and error space coherences are computed from the flux sensitivities of the dressed cavity transition frequencies as functions of the transmon-drive detuning, Eq. (4.93); at the common minimum the logical state |+L> is then propagated for 300 microseconds along 1/f flux noise tra- jectories under the Lindblad master equation, with and without the drive, and the decoherence-only infidelity of Eq. (4.94) and the Wigner function of the cavity state are reported. Calculation time: minutes

## Physical / mathematical content

- Quantum-technology examples. The files in this area model cavity QED, transmon qubits, NV centres, and related effective Hamiltonians. The recurring mathematics is finite-dimensional quantum dynamics with ladder operators, rotating-wave-style couplings, anharmonic oscillator terms, avoided crossings, and coherent control in coupled few-mode systems.
- The effective hardware model is a weakly anharmonic oscillator. Duffing nonlinearity breaks equal level spacing and allows qubit-like addressability within a truncated bosonic ladder.

## Numerical / algorithmic content

- An eigenvalue problem is solved or analysed, so the file is extracting spectra, stationary states, avoided crossings, or modal structure from the effective Hamiltonian or superoperator.
- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The Lindblad propagation is averaged over 100 windowed 1/f flux-noise trajectories from `pink_noise`, with the propagators tabulated on a 201-point grid of the transmon frequency offset and looked up per time step; the Wigner functions of the final cavity states are evaluated point by point on a 71x71 phase-space grid with `wigner_fock`, and the unit integral of the initial state on that grid is checked.
- The file also defines local helper function(s): `dressed_ens()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Implementation structure

- Binomial bosonic code |0L>=(|0>+|4>)/sqrt(2), |1L>=|2> in a cavity
- dispersively coupled to a flux-tunable transmon ancilla, and the
- protection of its coherences from 1/f flux noise by a Stark-assis-
- ted flux-noise evasion (SAFE) drive on the transmon, Sec. 4.4.1 and
- Fig. 4.4(a,b) of Yunwei Lu's PhD thesis (Northwestern University,
- 2026). The flux noise dephasing rates of the code and error space
- coherences are computed from the flux sensitivities of the dressed
- cavity transition frequencies as functions of the transmon-drive
- detuning, Eq. (4.93); at the common minimum the logical state |+L>
- is then propagated for 300 microseconds along 1/f flux noise tra-
- jectories under the Lindblad master equation, with and without the
- drive, and the decoherence-only infidelity of Eq. (4.94) and the

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `hamiltonian()`, `assume()`, `operator()`, `state()`, `int2str()`, `dressed_ens()`, `median()`, `num2str()`, `any()`, `kfigure()`, `subplot()`, `semilogy()`, `kxlabel()`.
