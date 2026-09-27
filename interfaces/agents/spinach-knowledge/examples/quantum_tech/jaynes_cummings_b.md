# examples/quantum_tech/jaynes_cummings_b.m

- Signature: `jaynes_cummings_b()`

## Purpose

A time-domain Jaynes–Cummings simulation of a spin coupled to an electromagnetic cavity mode. The initial state combines transverse spin magnetisation with an empty cavity; the trajectories report the spin's `Lx` and the cavity-field quadrature, and the example also plots cavity-level populations. Calculation time: seconds.

## Physical / mathematical content

- The system uses an electron spin and a five-level cavity mode, resonant at the electron frequency, with exchange coupling between them. The dynamics illustrate excitation exchange in the Jaynes–Cummings model.

## Numerical / algorithmic content

- Spinach constructs the system in the `sphten-liouv` formalism with no basis approximation and propagates it through the cavity device context. The trajectory is sampled at 251 points over 2.5 μs.

## Implementation structure

- The source sets `sys.magnet=0.33`, uses isotopes `{'E','C5'}`, and sets the cavity exchange to `2.828e6`. The initial state is `state(...,{'Lx','BL1'},{1,2}) + state(...,{'E','BL1'},{1,2})/2`; the sequence uses spin `E`, offset `5e6`, sweep `1e8`, and `npoints=251`. It projects the trajectory onto the spin `Lx`, the cavity quadrature `(C-A)/2i`, and the `BL1`–`BL3` cavity-level populations.
