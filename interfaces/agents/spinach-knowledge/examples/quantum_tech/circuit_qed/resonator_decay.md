# examples/quantum_tech/circuit_qed/resonator_decay.m

- Signature: `resonator_decay()`
- Source: [resonator_decay.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/quantum_tech/circuit_qed/resonator_decay.m)

## Purpose and open-system model

This example propagates a thermally damped microwave resonator and compares photon-number decay from a Fock state and a coherent state. In the finite-temperature bosonic relaxation model, emission and thermal absorption drive the mode toward its Bose–Einstein occupation. The finite Hilbert-space cutoff can distort that weak thermal channel, as the source comments explicitly caution.

The mode is truncated to five Fock levels (`C5`), with frequency 6.02 GHz, energy-relaxation lifetime 10 ns, and bath temperature 0.050 K. The code computes the thermal occupation as `n_eq = 1 / (exp(h f / (k_B T)) - 1)`, using the SI values of Planck’s and Boltzmann’s constants. It specifies a 5 μs pure-dephasing time and derives the total coherence time including thermal damping. The frequency/lifetime inputs and temperature are interpreted by the Spinach mode model; frequency-like inputs are supplied in Hz, lifetimes and times in seconds, and temperature in kelvin.

## Initial states, propagation, and observables

The initial conditions are the highest retained Fock level (`BL5`, four photons in this five-level basis) and a coherent state with amplitude 1.5. The harmonic-mode Hamiltonian is combined with Spinach’s finite-temperature bosonic relaxation superoperator; the Liouville-space generator is `G = -iH + R`. The script propagates both density operators on a 1 ns grid for 100 steps (0–100 ns). It records the populations of all five levels and obtains mean photon number by weighting those populations with 0 through 4.

For energy-damping rate `κ = 1/T1`, the comparison curve is `n(t) = (n(0) - n_eq) exp(-κ t) + n_eq`. The source checks both mean-photon trajectories against this curve to 5×10^-3, checks that the transient population maxima descend in Fock-level order, and requires the final Fock-run ground-state population to be within 2×10^-3 of `1/(1+n_eq)`. These are in-code numerical checks, not independently rerun results. The plotted population traces are simulated trajectories, not measured resonator data.

The source header attributes its model and parameters to the matching example in the paraqeet package.

## Scope

This is a five-level thermal open-system model. Its late-time state and accuracy are subject to that truncation; the script does not establish performance for an untruncated resonator or report an experimental decay measurement.
