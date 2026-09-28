# examples/nmr_liquids/noesy_sucrose.m

- Signature: `noesy_sucrose()`

## Purpose

NOESY spectrum of sucrose (magnetic parameters computed with DFT). Calculation time: minutes.

## Physical / mathematical content

- Sucrose magnetic parameters are read from a vacuum-DFT log with a minimum retained J-coupling of 1.0 Hz. The proton spin system is simulated at 5.9 T.
- Redfield relaxation uses the IME equilibrium option, 298 K, a 200 ps correlation time, and kite retention; the NOESY mixing time is 0.5 s.

## Numerical / algorithmic content

- The spherical-tensor Liouville basis uses IK-2 with scalar-coupling connectivity and proximity level 3. Greedy handling is enabled, Krylov propagation is disabled, and the proximity cutoff is 4.0.
- The 1H acquisition uses an 800 Hz offset, 1700 Hz sweep in both dimensions, 512 acquired points and 2048-point zero filling per dimension. Square-cosine apodisation and States processing precede the two Fourier transforms.

## Implementation structure

- Parse the sucrose DFT log, build the proton basis and configure the NOESY sequence.
- Simulate, apodise cosine and sine FIDs, combine the States signal, Fourier transform F2 and F1, then plot the negative real spectrum.
