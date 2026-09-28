# examples/parahydrogen/pasadena_propanal.m

- Signature: `pasadena_propanal()`

## Purpose

Simulates a PASADENA spectrum for parahydrogenation of acrolein to propanal. The source gives a calculation time of seconds.

## Physical / mathematical content

The six-spin system models the propanal product at 7.05 T using the specified proton chemical shifts and scalar couplings. It uses the spherical-tensor Liouville formalism with no basis approximation and the `S3` and `S2` symmetry groups on spins 1–3 and 4–5. The initial state is the `Lz` product on spins 1 and 4.

## Numerical / algorithmic content

The proton liquid-state acquisition uses a `pi/4` pulse, 500 ppm offset, 1000 ppm sweep, 1024 points, and 8192-point zero filling. The FID is Gaussian-apodised with parameter 10, Fourier transformed, and plotted.

## Implementation structure

The script defines six proton isotopes, the field and scalar interaction data, creates the symmetry-adapted basis, and runs `liquid` with `hp_acquire` before signal processing.
