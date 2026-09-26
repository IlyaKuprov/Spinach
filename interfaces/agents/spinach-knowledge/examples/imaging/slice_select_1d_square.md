# examples/imaging/slice_select_1d_square.m

- Signature: `slice_select_1d_square()`

## Purpose

Slice selection example using a one-dimensional sample and a square slice selection pulse in the presence of diffusion and flow. Calculation time: seconds. Ahmed Allami Ilya Kuprov

## Physical / mathematical content

- MRI and spectroscopic-imaging examples. These files combine gradient terms, spatial encoding, diffusion, slice selection, k-space sampling, and Fourier reconstruction, generally within Fokker-Planck or explicit spatial-grid descriptions.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- Slice selection example using a one-dimensional sample and
- a square slice selection pulse in the presence of diffusion
- and flow.
- Calculation time: seconds.
- Ahmed Allami
- Ilya Kuprov
- Isotopes
- Magnetic induction
- Chemical shifts
- Relaxation model
- Basis set
- Spinach housekeeping
