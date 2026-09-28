# examples/quantum_tech/transmon_rabi_leakage.m

- Signature: `transmon_rabi_leakage()`

## Purpose

Rabi dynamics of a driven four-level transmon in the Duffing approximation, including leakage into the second and third excited states. The resonant drive is part of the rotating-frame Hamiltonian, and all four level populations come from a single trajectory. Calculation time: seconds.

## Model and parameters

The model is a T4 transmon at zero rotating-frame frequency with anharmonicity -250 MHz, in the Zeeman-Hilbert formalism without approximation. A resonant 25 MHz drive is included in the Hamiltonian, and the initial state is BL1.

## Calculation

The code propagates one trajectory with a 1 ns step for 400 ns, evaluates populations in BL1 through BL4 at every point, and plots all four traces. Population outside the lowest two levels displays leakage into the higher transmon states.
