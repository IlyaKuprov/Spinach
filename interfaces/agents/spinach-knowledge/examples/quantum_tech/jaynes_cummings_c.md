# examples/quantum_tech/jaynes_cummings_c.m

- Signature: `jaynes_cummings_c()`

## Purpose

A time-domain Jaynes–Cummings simulation of two exchange-coupled electrons, each coupled to the same electromagnetic cavity mode. The initial state contains transverse spin magnetisation and an empty cavity; the calculation detects the summed spin `Lx` signal and the cavity-field quadrature. Calculation time: seconds.

## Physical / mathematical content

- The cavity is resonant with the electrons, and the two spins have distinct cavity exchange couplings. Their mutual scalar exchange coupling is included, so the spin and cavity dynamics evolve together.

## Numerical / algorithmic content

- The source builds the system in Spinach's `sphten-liouv` formalism with no basis approximation, then propagates it through the cavity device context. The trajectory uses 251 points over 2.5 μs.

## Implementation structure

- Parameters include `sys.magnet=0.33`, isotopes `{'E','E','C5'}`, electron–electron scalar coupling `5e6`, cavity exchange couplings `2.828e6` and `2.728e6`, sequence offset `5e6`, sweep `1e8`, and `npoints=251`. Both spins start with `Lx` contributions paired with cavity state `BL1`; the detected spin operator sums their `Lx` operators, and the cavity signal uses `(C-A)/2i` on mode 3.
