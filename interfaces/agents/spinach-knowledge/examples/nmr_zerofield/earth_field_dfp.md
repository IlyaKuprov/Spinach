# examples/nmr_zerofield/earth_field_dfp.m

- Signature: `earth_field_dfp()`

## Purpose

Simulates Earth-field NMR of 2,6-difluoropyridine, reproducing the simulated spectra in Figure 7 of [the cited paper](https://doi.org/10.1016/j.jmr.2023.107540) without the weighted addition of the uncoupled 1H signal. The source estimates seconds. Authors: Adam Altenhof (adamaltenhof@gmail.com) and Derrick Kaseman (kaseman1@llnl.gov).

## Spin system and sequence

The field is set from 2312.45 Hz for 1H. The model has labeled H3, H4, H5, F2, F6, and N1 spins (1H, 19F, and 14N), with shifts 6.98, 8.06, 6.98, -70.69, -70.69, and 0. The experimental coupling values include F-N couplings of 37.35 Hz, H-N couplings set to (1.03+0.71)/2 Hz, H-F couplings of -2.47, 8.08, and 1.29 Hz, an F-F coupling of -12.23 Hz, H-H couplings of 7.92 Hz, and an H3-H5 coupling of 0.55 Hz. The source notes that Table 2 is unclear for the H-N values.

The untruncated sphten-liouv basis uses T1/T2 relaxation, diagonal retention, and zero equilibrium. Sequence settings are sweep=1e4 Hz, 60000 points, zero-fill to 2^19, zero offset, 1H detection in Hz, uniaxial detection, and flip angle pi/2. For the 14N relaxation sweep, its T1 and T2 rates take values 0, 10, 100, 500, 1e3, 1e4, 1e5, and 1e6 Hz; the other five spins use 0.22 Hz. The source comments that a GPU is needed, while the GPU-enable line is commented out.

## Simulation and processing

For each 14N rate, the example calls the liquid simulator in the lab frame, Fourier transforms and normalizes the real spectrum, then plots the 2303-2322 Hz region.
