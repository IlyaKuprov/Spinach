# examples/nmr_solids/static_powder_nqi_a.m

- Signature: `static_powder_nqi_a()`

## Purpose

Simulates the static 14N powder pattern of L-valyl-L-alanine, using the large orientation grid specified in the example to reproduce Figure 5 of O'Dell and Ratcliffe. The source estimates a calculation time of minutes. Reference: [O'Dell and Ratcliffe](https://doi.org/10.1016/j.cplett.2011.08.030).

## Spin system and interactions

The model contains two 14N spins at 21.1 T. Their quadrupolar interaction matrices are generated with `eeqq2nqi`, using coupling magnitudes 1.24 and 3.06 MHz and asymmetries 0.22 and 0.40, respectively. The basis is the full Zeeman Hilbert-space basis (`zeeman-hilb`, no approximation).

## Simulation and processing

The acquisition uses 14N, a 6 MHz sweep, 512 points, a 2048-point zero-fill, and the `icos_2ang_163842pts` powder grid. The frequency axis is in MHz and inverted. Initial and detection states are both the 14N `L+` state. The powder FID is apodised with an exponential parameter of 6, Fourier transformed, and plotted.
