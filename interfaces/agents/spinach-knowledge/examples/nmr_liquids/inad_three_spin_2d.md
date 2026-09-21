# examples/nmr_liquids/inad_three_spin_2d.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_liquids/inad_three_spin_2d.m`
- Signature: `inad_three_spin_2d()`
- Total lines: 81

## Purpose

Example of a 2D-INADEQUATE spectrum of a generic three-spin system. Calculation time: seconds. Theresa Hune Christian Griesinger

## Physical / mathematical content

- Liquid-state NMR examples. The physics is scalar-coupling-mediated coherence transfer in weakly or moderately coupled spin systems, often in Liouville space. Typical mechanisms include INEPT-style polarisation transfer, J-refocusing, phase cycling, indirect evolution, and multidimensional detection.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.
- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Implementation structure

- Example of a 2D-INADEQUATE spectrum of a
- generic three-spin system.
- Calculation time: seconds.
- Theresa Hune
- Christian Griesinger
- Magnetic field (700 MHz)
- Generic three-spin system
- Formalism and basis set
- Spinach housekeeping
- Generate isotopomers
- Sequence parameters
- Preallocate the answer

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `dilute()`, `liquid()`, `apodisation()`, `fftshift()`, `kfigure()`, `scale_figure()`, `plot_2d()`, `kylabel()`.
