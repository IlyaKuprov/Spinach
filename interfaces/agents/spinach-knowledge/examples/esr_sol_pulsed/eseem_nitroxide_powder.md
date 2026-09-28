# examples/esr_sol_pulsed/eseem_nitroxide_powder.m

- Signature: `eseem_nitroxide_powder()`

## Purpose

Powder-averaged two-pulse ESEEM on a 14N nitroxide radical. Time-domain simulation in Liouville space with powder averaging over a finite grid. Set to reproduce Figure 4a in http://dx.doi.org/10.1063/1.453532; ideal pulses are assumed. Calculation time: seconds.

## Physical / mathematical content

- The system contains 14N and an electron at 0.3249 T. Nitrogen self-coupling eigenvalues are [-0.4, -1.6, 2.0]×10^5 and electron–nitrogen coupling eigenvalues are [2, 2, 2]×10^6; both interactions have zero Euler angles.
- Powder averaging uses `rep_2ang_400pts_sph`.

## Numerical / algorithmic content

- The sequence uses 512 points at a 200 ns timestep and 2048-point zero filling. The signal is mean-subtracted and exponentially apodised with parameter 5 before Fourier transformation and `fftshift`.
- The frequency axis uses an interpulse-delay increment of half the timestep.

## Implementation structure

- Create the spin system in the `sphten-liouv` basis without approximation, with trajectory-level SSR disabled; call `powder` with `@eseem` in the `esr` context.
- Plot the real apodised time-domain signal and real spectrum.
