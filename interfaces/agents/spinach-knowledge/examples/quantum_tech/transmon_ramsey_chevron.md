# examples/quantum_tech/transmon_ramsey_chevron.m

Source: [examples/quantum_tech/transmon_ramsey_chevron.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/quantum_tech/transmon_ramsey_chevron.m)

## What it models

A closed, coherently controlled three-level transmon in the Duffing approximation. The plotted chevron is a simulated Ramsey sequence over detuning and free-evolution time, not an experimental measurement; the source estimates calculation time in seconds.

## Hamiltonian and parameters

The source sets the field to zero, uses a T3 mode with rotating-frame frequency 0 and anharmonicity -260e6 Hz, and selects the Zeeman-Hilbert formalism without basis approximation. The cavity/Duffing Hamiltonian supplies the anharmonic drift H0. For each detuning Δ in a 256-point grid from -20e6 to 20e6 Hz, the evolution Hamiltonian is H0 + 2*pi*Δ*N, symmetrised in the source; N is the transmon number operator. The 256 free-evolution times span 0 to 1.0e-6 s.

## Pulse sequence and detection

The initial BL1 state is acted on by a nominal pi/2 propagator made from the transmon quadrature (C+A)/2. At every detuning and time sample, the code propagates under the detuned drift, applies the same pi/2 propagator again, and detects BL2 population. The image uses time in seconds on the horizontal axis and detuning in Hz on the vertical axis, with colour encoding the calculated population.

No dissipative relaxation or dephasing term, rotational diffusion, correlation spectrum, cross-correlation, or secular approximation is specified. The sequence therefore illustrates ideal coherent Ramsey fringes; it does not establish experimental agreement or a coherence-time benchmark.
