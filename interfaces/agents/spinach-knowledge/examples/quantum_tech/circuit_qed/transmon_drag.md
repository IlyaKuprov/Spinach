# examples/quantum_tech/circuit_qed/transmon_drag.m

- Signature: `transmon_drag()`
- Source: [transmon_drag.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/quantum_tech/circuit_qed/transmon_drag.m)

## Purpose and physical model

The example compares a plain resonant Gaussian pulse, its analytic DRAG correction, and a numerically tuned DRAG pulse on a three-level Duffing transmon. Because the computational transition is anharmonic rather than isolated from the next level, a short resonant pulse can leak population into the second excited state. DRAG adds an envelope-derivative contribution to the orthogonal IQ quadrature, scaled by a detuning parameter, to reduce that leakage while implementing a qubit-subspace π/2 (90°) rotation.

The transmon has a Duffing-oscillator drift with a 4.8 GHz 0–1 frequency and −200 MHz self-anharmonicity. These mode inputs are in Hz and are converted internally to angular-frequency units. The driven Hamiltonian adds the transmon X quadrature operator (`C + A`) multiplied by the carrier-modulated Gaussian plus its DRAG derivative quadrature. The pulse lasts 20 ns, with Gaussian width `sigma = t_dur/8 = 2.5 ns`, and is integrated on 2,000 midpoint steps of 10 ps. The analytic detuning is `2 × 2π × anharm`, i.e. approximately −2.5133×10^9 rad/s for the coded anharmonicity.

## Pulse settings and computed quantities

The three rows in `pulse_params` are the plain Gaussian, analytic DRAG, and optimised DRAG cases; columns are amplitude, detuning, local-oscillator angular frequency, phase, and DRAG-quadrature switch. The first two use amplitude 3.0×10^8 rad/s and `2π × 4.8 GHz`, with zero phase; only the analytic row enables the derivative quadrature. The optimised row uses amplitude 2.455066180403117×10^8 rad/s, detuning −2.492899680375753×10^9 rad/s, oscillator frequency 3.015929426382811×10^10 rad/s, phase −1.886372179487061×10^-4 rad, and enables DRAG.

For each pulse, the ket starts in the ground state and is propagated with the laboratory-frame Hamiltonian. The resulting propagator is transformed into the drift frame before comparison with the target π/2 rotation. The reported score is the three-level process-fidelity expression `abs(mean(diag(target' * U)))^2`; it includes the leakage level in the comparison. The source’s embedded reference vector is `[0.969588, 0.973797, 0.994239]`, with a 0.001 tolerance. It also checks that the final population in level 3 for analytic DRAG is less than half the plain-pulse value. These are simulated reference/check values encoded in the MATLAB example, not experimental device measurements or results rerun for this draft.

The source header attributes its model and parameters to the matching example in the paraqeet package.

## Scope

This is closed-system three-level pulse simulation. It illustrates how the derivative quadrature changes leakage and the computed gate score for the coded parameters; it does not establish a hardware gate fidelity or a globally optimal DRAG pulse.
