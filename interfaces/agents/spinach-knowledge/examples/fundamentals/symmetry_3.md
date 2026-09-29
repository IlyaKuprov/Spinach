# examples/fundamentals/symmetry_3.m

- MATLAB implementation: [examples/fundamentals/symmetry_3.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/fundamentals/symmetry_3.m)

- Signature: `symmetry_3()`
- Source: [examples/fundamentals/symmetry_3.m](../../../../../examples/fundamentals/symmetry_3.m)

## Purpose

Produces a simulated 1H NMR spectrum for the source's eight-proton valine model using the fully symmetric irreducible representation of S3(x)S3.

## Spin model and basis

The field is 11.7 T and all eight spins are protons. Zeeman scalar values in spin order are 3.5950, 2.2580, 1.0270 (three entries), and 0.9760 (three entries). The scalar couplings are 4.34 for spins 1–2 and 7.00 between spin 2 and each of spins 3–8; the source also assigns 0.00 to entry (8,8). Two S3 groups act on spins [3 4 5] and [6 7 8]. The basis is spherical-tensor Liouville with approximation `none`.

## Acquisition and processing

Both the initial state and coil are `state(spin_system,'L+','1H')`, and the decoupling list is empty. Liquid-state acquisition uses 8192 points, sweep 2500 Hz, and offset 1000 Hz. The FID receives exponential apodisation `{'exp',5}`, is zero-filled to 65536 points, then Fourier-transformed and plotted using the real spectrum in ppm with the axis inverted.

## Scope

The page records the model and processing parameters present in the source. It does not claim observed peak positions, a measured spectrum, or a run result.
