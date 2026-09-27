# examples/nmr_liquids/noesy_strychnine.m

- Signature: `noesy_strychnine()`

## Purpose

NOESY spectrum of strychnine. Calculation time: minutes.

## Physical / mathematical content

- The proton spin system is strychnine. Redfield relaxation uses the IME equilibrium option, 298 K, a 200 ps correlation time, and kite retention. The sequence has a 0.5 s mixing time.
- The simulated two-dimensional spectrum is derived from the NOESY cosine and sine FIDs using States combination and Fourier transformation in both dimensions.

## Numerical / algorithmic content

- The spherical-tensor Liouville basis uses IK-2 with scalar-coupling connectivity and proximity level 3. Greedy handling is enabled, Krylov propagation disabled, and the proximity cutoff is 4.0.
- The 1H spectrum uses a 1200 Hz offset, 2500 Hz sweep in each dimension, 512 acquired points and 2048-point zero filling per dimension. Square-cosine apodisation is applied to both States components.

## Implementation structure

- Create the strychnine proton spin system at 5.9 T, construct the basis and set the 0.5 s NOESY mixing period.
- Simulate the sequence, apodise the cosine and sine FIDs, combine the States signal, Fourier transform F2 and F1, and plot the negative real spectrum.
