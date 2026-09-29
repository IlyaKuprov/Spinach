# examples/relaxation_theory/scalar_relaxation_1.m

[Source file](https://github.com/IlyaKuprov/Spinach/blob/main/examples/relaxation_theory/scalar_relaxation_1.m) · Signature: `scalar_relaxation_1()`

## Purpose and physical setting

Builds a relaxation superoperator for scalar relaxation of the first kind in a two-proton model with a noisy J-coupling. The source relates the example to aziridines, where slow nitrogen inversion jitters scalar couplings on a millisecond time scale, and cites [the described effect](https://doi.org/10.1002/ange.201410271). The model itself contains the two listed proton spins; it is not a simulated nitrogen-inversion trajectory.

## Model and output

The source sets `sys.magnet=11.75` and `inter.zeeman.scalar={0.0 2.0}`; units for these scalar values are not identified. It uses the complete `sphten-liouv` basis without approximation and configures `inter.relaxation={'SRFK'}`, `inter.rlx_keep='kite'`, and zero equilibrium. The SRFK correlation-time input is `[1.0 1e-3]`, with off-diagonal modulation-depth entry `15.0`; the source does not label units for these values. After system and basis construction, it evaluates `relaxation(spin_system)` and displays the superoperator's nonzero pattern with `spy`. It does not set up an acquisition, pulse sequence, or spectrum calculation.
