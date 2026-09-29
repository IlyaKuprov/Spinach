# examples/quantum_tech/transmon_rabi_leakage.m

Source: [examples/quantum_tech/transmon_rabi_leakage.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/quantum_tech/transmon_rabi_leakage.m)

## What it models

A closed, coherently driven four-level transmon in the Duffing approximation. The four truncated ladder states are reported as BL1 through BL4, so the higher two populations show leakage beyond the computational pair in this model. This is a calculated trajectory, not experimental data; the source estimates calculation time in seconds.

## Hamiltonian and parameters

The source sets the field to zero and declares one T4 mode, zero rotating-frame frequency, and anharmonicity -250e6 Hz. It uses the Zeeman-Hilbert formalism with no basis approximation and builds the cavity/Duffing drift Hamiltonian. A resonant quadrature drive is added as 2*pi*25e6*(C+A)/2, where C and A are the transmon ladder operators; the 25 MHz source parameter is converted to angular frequency in the Hamiltonian.

## State, propagation, and plot

The initial state is BL1. Spinach propagates a single trajectory with 1 ns steps for 400 ns. At each point the code evaluates the BL1, BL2, BL3, and BL4 populations using the corresponding state operators and plots all four against time in ns. The figure visualises leakage during ideal coherent Rabi dynamics; the source includes no relaxation, decoherence, rotational diffusion, correlation spectrum, cross-correlations, or secular approximation, and reports no comparison with an experiment.
