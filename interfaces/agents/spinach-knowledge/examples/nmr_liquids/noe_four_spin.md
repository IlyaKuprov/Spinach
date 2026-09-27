# examples/nmr_liquids/noe_four_spin.m

- Signature: `noe_four_spin()`

## Purpose

Inversion-recovery NOE effect spectrum on a simple four-spin system, with the rightmost proton signal inverted and a pulse-acquire experiment performed after a very long (five seconds) mixing time. Sequential NOE hops with alternating signs are clearly visible in the result. Calculation time: seconds.

## Implementation

The four 1H spin system uses Redfield relaxation, the Di Bari equilibrium state, a 298 K temperature, and a 200 ps correlation time. The code prepares the state with the rightmost spin inverted, evolves its deviation from equilibrium for five seconds under the relaxation superoperator, then applies a 1H pulse-acquire experiment. It exponentially apodises and Fourier transforms the resulting FID for plotting.
