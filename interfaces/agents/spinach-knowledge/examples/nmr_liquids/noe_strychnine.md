# examples/nmr_liquids/noe_strychnine.m

- Signature: `noe_strychnine()`

## Purpose

Inversion-recovery NOE effect spectrum on strychnine, with the rightmost proton signal inverted and a pulse-acquire experiment performed after a 500 ms mixing time. Calculation time: minutes.

## Physical / mathematical content

- The strychnine proton spin system is treated with scalar-coupling connectivity and Redfield relaxation. The relaxation parameters are the Di Bari equilibrium convention, 298 K, and a 200 ps correlation time.
- The initial state is formed by inverting spin 9 relative to thermal equilibrium. After 500 ms of relaxation evolution, the unperturbed equilibrium state is subtracted before acquisition, isolating the NOE difference signal.

## Numerical / algorithmic content

- Spinach uses the spherical-tensor Liouville basis with the IK-2 approximation, scalar-coupling connectivity, proximity level 3, and a 4.0 proximity cutoff; Krylov propagation is disabled.
- A pulse-acquire calculation uses 8192 points, zero-filled to 65536, with a 6500 Hz sweep and 2800 Hz offset. The FID receives 6 Hz exponential apodisation before the Fourier transform.

## Implementation structure

- Build the strychnine spin system and Redfield relaxation superoperator, calculate equilibrium, and construct the spin-9-inverted state.
- Evolve for 0.5 s, subtract equilibrium, then run pulse-acquire with a proton coil and 90-degree `Ly` pulse. Invert the plotted frequency axis and display the real spectrum.
