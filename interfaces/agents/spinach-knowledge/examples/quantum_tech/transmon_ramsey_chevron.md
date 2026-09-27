# examples/quantum_tech/transmon_ramsey_chevron.m

- Signature: `transmon_ramsey_chevron()`

## Purpose

Ramsey chevron of a three-level transmon in the Duffing approximation. A nominal pi/2 pulse prepares a coherence, and detuning during free evolution produces Ramsey fringes. Calculation time: seconds.

## Model and parameters

The T3 transmon is in the rotating frame with anharmonicity -260 MHz and uses the Zeeman-Hilbert formalism without approximation. The detuning grid contains 256 points from -20 to 20 MHz; free-evolution times contain 256 points from 0 to 1 microsecond.

## Calculation

A nominal pi/2 propagator prepares the initial BL1 state. For each detuning and time, the code propagates under the anharmonic Hamiltonian plus the number-operator offset, applies the same final pi/2 pulse, and detects BL2 population. The result is plotted as a detuning-versus-time Ramsey chevron.
