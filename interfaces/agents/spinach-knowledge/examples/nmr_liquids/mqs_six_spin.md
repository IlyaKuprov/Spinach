# examples/nmr_liquids/mqs_six_spin.m

- Signature: `mqs_six_spin()`
- Source: [`examples/nmr_liquids/mqs_six_spin.m`](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_liquids/mqs_six_spin.m)

## Purpose and spin system

A simulated multiple-quantum NMR experiment on six coupled 1H spins; the source estimates minutes of calculation time. The field parameter is 14.1, and the six chemical-shift entries are `[-1.0, -0.5, 0.0, 0.3, 0.7, 1.1]`. The scalar-coupling assignments are: 3J, pairs (1,2)=8.77, (2,3)=7.40, (3,4)=7.40, (4,5)=8.77, (5,6)=8.77; 4J, pairs (1,3), (2,4), (3,5), and (1,6), each 2.60; 5J, pairs (1,4), (2,5), (2,6), and (4,6), each 2.00; 6J, pairs (1,5) and (3,6), each 2.00. The file also sets `inter.coupling.scalar{6,6}=0`. The source does not annotate units for the shifts or coupling entries. The basis is the full `sphten-liouv` basis without approximation.

## Coherence experiment and spectrum

The initial proton state is `Lz`, the detected state is `L+`, and the pulse angle is `pi/2`. The selected coherence orders are `[+6 -1]`; the sequence uses two 1H dimensions, zero offsets (explicitly in Hz), sweeps `[2400 2400]` Hz, 512 points per dimension, and zero-filling to 2048 per dimension. The axis units are ppm. The two delays are each 0.1 s.

The example calls `liquid(spin_system,@mqs_refocus,parameters,'nmr')`, applies squared-cosine apodisation in both dimensions, and computes a zero-filled 2D FFT. It plots the absolute spectrum with axes labelled 1Q / ppm and 6Q / ppm. No relaxation model, experimental dataset, measured transfer rate, or DOI is specified in this source.
