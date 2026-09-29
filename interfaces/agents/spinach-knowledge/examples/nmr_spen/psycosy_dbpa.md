# examples/nmr_spen/psycosy_dbpa.m

- MATLAB implementation: [examples/nmr_spen/psycosy_dbpa.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_spen/psycosy_dbpa.m)

- Signature: `psycosy_dbpa()`
- Source: [`examples/nmr_spen/psycosy_dbpa.m`](../../../../../examples/nmr_spen/psycosy_dbpa.m)

## Spin model and spatial encoding

This PSYCOSY example represents the dibromopropionic-acid ring with four protons at 14.1 T. The source assigns shift entries [4.49, 3.9, 3.7, 4.2] and nonzero scalar-coupling entries (1,2)=11.3, (2,3)=10.1, and (1,3)=4.3; entries involving site 4 are explicitly set to zero. The model uses the spherical-tensor Liouville formalism, IK-2 approximation, proximal level 1, scalar-coupling connectivity, and greedy basis construction.

The `@psycosy` sequence runs through `imaging` on a 15 mm sample represented by 100 points, with derivative option `{'period',3}`. Diffusion and flow are zero. A zero-valued relaxation phantom and relaxation superoperator are supplied; the initial and receive states are `Lz` and `L+`.

## Sequence and processing

The sequence assignments are offset 2460, sweep 600, acquisition sizes [512, 512], zero-fill sizes [1024, 1024], and `axis_units='ppm'`. The mixing time is 25e-3 s and gradient amplitude 1e-2 T/m. The saltire chirp uses a 20-degree flip angle, 0.015 s duration, 0.05 s chirp-gradient duration, 10000 Hz sweep width, 250 pulse points, and smoothing factor 20. The source does not annotate units for the offset or the `sweep` assignment.

The two-dimensional FID is apodised with square-sine windows in both dimensions, transformed with a zero-filled 2D FFT, and plotted as its magnitude. The source header estimates minutes on an NVIDIA Tesla A100 and says CPU execution takes much longer; this is the source's estimate, not a run result.
