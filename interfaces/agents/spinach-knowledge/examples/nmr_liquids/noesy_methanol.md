# examples/nmr_liquids/noesy_methanol.m

- Signature: `noesy_methanol()`

## Purpose

NOESY spectrum of 13C methanol. J-couplings from Pecul and Helgaker, CSA tensors from DFT. Note the presence of cross-peaks between 13C doublet components. Calculation time: seconds.

## Physical / mathematical content

- The magnetic parameters are parsed from the methanol vacuum-DFT log; the OH proton is then removed, leaving the carbon and three methyl protons. The isotropic shifts of the four retained spins are set on resonance.
- Scalar couplings are assigned as 141 Hz between carbon and each proton, and −11 Hz between each pair of protons. Redfield relaxation uses the IME equilibrium option, 298 K, a 50 ps correlation time, and kite retention.
- The NOESY spectrum exhibits cross-peaks between the two components of the carbon-13 doublet.

## Numerical / algorithmic content

- A spherical-tensor Liouville basis without approximation is used, with greedy state-space handling and a 4.0 proximity cutoff.
- The NOESY mixing time is 0.5 s. The two dimensions use 300 Hz sweeps, 256 acquired points and 1024-point zero filling; both States components receive square-cosine apodisation before sequential Fourier transforms.

## Implementation structure

- Parse the methanol DFT output, remove the OH spin and its interaction entries, assign the field (14.1 T), chemical shifts, and listed couplings, then build the basis.
- Run the NOESY sequence, apodise cosine and sine FIDs, combine them as a States signal, Fourier transform both dimensions, and plot the negative real spectrum.
