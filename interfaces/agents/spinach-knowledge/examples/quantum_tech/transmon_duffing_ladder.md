# examples/quantum_tech/transmon_duffing_ladder.m

- Signature: `transmon_duffing_ladder()`

## Purpose

Duffing-model energy ladder of a weakly anharmonic transmon, showing how the transition frequencies separate as anharmonicity increases. Inspired by the transmon model of Koch et al., Phys. Rev. A 76, 042319 (2007). Calculation time: seconds.

## Model and parameters

A single T5 mode has frequency 5.0 GHz. The code forms its harmonic Hamiltonian in the lab-frame context, then adds the Duffing term using the CCAA operator. Anharmonicity is swept across 80 values from -400 to -50 MHz.

## Calculation

For each anharmonicity, the Hamiltonian eigenvalues are sorted and the four adjacent transition frequencies are computed. The plot shows those transitions against `-alpha/2pi` in MHz, with frequencies in GHz. The cited transmon model is Koch et al., Phys. Rev. A 76, 042319 (2007).
