# examples/imaging/slice_select_1d_shaped.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/imaging/slice_select_1d_shaped.m`
- Signature: `slice_select_1d_shaped()`
- Total lines: 88

## Purpose

Slice selection example using a one-dimensional sample and a shaped slice selection pulse in the presence of diffusion and flow. Calculation time: seconds. Ahmed Allami Ilya Kuprov

## Physical / mathematical content

- MRI and spectroscopic-imaging examples. These files combine gradient terms, spatial encoding, diffusion, slice selection, k-space sampling, and Fourier reconstruction, generally within Fokker-Planck or explicit spatial-grid descriptions.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- Slice selection example using a one-dimensional sample and
- a shaped slice selection pulse in the presence of diffusion
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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `pulse_shape()`, `relaxation()`, `state()`, `imaging()`, `apodisation()`, `fftshift()`, `ifftshift()`, `kfigure()`, `plot_1d()`.
