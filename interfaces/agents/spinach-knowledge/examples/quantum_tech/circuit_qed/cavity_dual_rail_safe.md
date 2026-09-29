# examples/quantum_tech/circuit_qed/cavity_dual_rail_safe.m

Source: https://github.com/IlyaKuprov/Spinach/blob/main/examples/quantum_tech/circuit_qed/cavity_dual_rail_safe.m

## Purpose and model

This example asks whether one Stark-assisted flux-noise evasion (SAFE) drive can reduce dephasing of two dual-rail qubits sharing a flux-tunable transmon ancilla, as in Section 4.4.2 and Figure 4.4(c) of Yunwei Lu's 2026 Northwestern University PhD thesis. The modeled phase-sensitive transitions are the single-photon states of the two transmon-coupled cavity rails relative to their ground state. The finite model contains a three-level transmon (`T3`) and two three-level cavity modes (`C3`, `C3`); it does not include explicit reference rails or a full hardware device.

The transmon anharmonicity is −200 MHz. The two rail couplings are 50 and 60 MHz and their detunings are 1.414 and 2.0 GHz, respectively. The effective drift has the transmon self-anharmonicity and a transmon-number/cavity-number cross shift for each rail. The off-resonant SAFE drive Stark-dresses the transmon to change the dressed single-photon transitions’ flux slopes. Its flux-noise operator weights transmon number, both rail numbers, and the two transmon–rail number products by their sensitivities. The dispersive shifts and sensitivities are derived from those parameters. The SAFE drive amplitude is 10 MHz; its detuning is scanned from −20 to −50 MHz in 0.5 MHz steps. The flux-noise amplitude is `1e-5` flux quanta, and the transmon frequency sensitivity is `2*pi*6e9` rad/s per flux quantum.

## Dephasing objective and interpretation

At each detuning, the code diagonalises the rotating-frame Hamiltonian (thesis Eq. 4.16) with symmetric frequency perturbations of `±2*pi*1e3` rad/s and obtains dressed-energy sensitivities for the ground and single-photon states. It converts each differential sensitivity into a flux-noise dephasing rate using the thesis expressions in Eqs. (4.93) and (4.95), including the source's fixed logarithmic factor `|ln(omega_ir*t)| = 4`. The driven rates are compared with the undriven dispersive estimates. A common operating point is selected from the two minima; the script checks that neither minimum lies at the scan boundary, that the minima are within 5 MHz of each other, and that the drive suppresses each rate by at least a factor of five.

This is a rate-versus-detuning calculation from an effective Hamiltonian, not a time-domain noise-trajectory simulation: the source defines no relaxation times, noise trajectories, or experimental observations. The plotted rates and suppression factors are model predictions for the stated parameter set. The source implements the thesis calculation but does not demonstrate experimental protection or convergence beyond its three-level truncations and detuning grid. No MATLAB output values are asserted here.
