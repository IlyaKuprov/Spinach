# examples/nmr_liquids/hsqc_cyprinol.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_liquids/hsqc_cyprinol.m`
- Signature: `hsqc_cyprinol()`
- Total lines: 80

## Purpose

HSQC spectrum of cyprinol with natural abundance of 13C isotope. Calculation time: seconds

## Physical / mathematical content

- Liquid-state NMR examples. The physics is scalar-coupling-mediated coherence transfer in weakly or moderately coupled spin systems, often in Liouville space. Typical mechanisms include INEPT-style polarisation transfer, J-refocusing, phase cycling, indirect evolution, and multidimensional detection.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.
- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Implementation structure

- HSQC spectrum of cyprinol with natural abundance of 13C isotope.
- Calculation time: seconds
- Spin system
- Magnet field
- Algorithmic options
- Basis set
- Sequence parameters
- Create the spin system structure
- Generate isotopomers
- Preallocate the answer
- Loop over isotopomers
- Build the basis

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `cyprinol()`, `create()`, `dilute()`, `basis()`, `liquid()`, `apodisation()`, `fftshift()`, `conj()`, `kfigure()`, `scale_figure()`, `plot_2d()`.
