# examples/parahydrogen/pasadena_ethylbenzene.m

- Signature: `pasadena_ethylbenzene()`

## Purpose

Simulates a PASADENA spectrum for parahydrogenation of styrene to ethylbenzene, with the stated aim of reproducing the top trace of Figure 5 ([doi:10.1039/b914188j](https://doi.org/10.1039/b914188j)). The source gives a calculation time of seconds.

## Physical / mathematical content

The model contains ten proton spins at 7.05 T, with chemical shifts and scalar couplings assigned to the product. It uses the spherical-tensor Liouville formalism, an IK-2 basis approximation connected through scalar couplings, and permutation symmetries `S3` for spins 1–3 and `S2` for spins 4–5. The initial state is the `Lz` product on spins 1 and 4.

## Numerical / algorithmic content

A liquid-state acquisition simulates a proton FID with a `-pi/4` pulse, 500 ppm offset, 1000 ppm sweep, 1024 points, and 8192-point zero filling. Gaussian apodisation (parameter 10) precedes the Fourier transform and spectrum plot.

## Implementation structure

The script defines the ten isotopes, field, shifts, and scalar couplings; constructs the symmetry-reduced basis; sets the proton acquisition parameters; and calls `liquid` with `hp_acquire` before processing the FID.
