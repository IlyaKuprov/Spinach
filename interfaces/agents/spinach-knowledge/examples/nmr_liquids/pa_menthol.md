# examples/nmr_liquids/pa_menthol.m

- Signature: `pa_menthol()`
- Source: [MATLAB example](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_liquids/pa_menthol.m)

## Purpose

Simulates a liquid-state menthol 1H NMR FID and illustrates effects labelled as bad Z1 and Z2 magnet shims. It is a pulse-acquire calculation using `liquid(...,@acquire,...,'nmr')`, not an INADEQUATE, inversion-recovery, NOE, or NOESY sequence.

## System and basis

The example loads `sys` and `inter` from `menthol.mat`; their numerical shifts, couplings, isotope inventory, and field are not stated in this MATLAB file. It builds a scalar-coupling-connected `sphten-liouv` basis with `IK-2` approximation, proximal level 1, projection +1, and three S3 symmetry groups. The source enables the greedy algorithm. Zero track elimination is also explicitly enabled with `zte` in `sys.enable`.

## Acquisition and processing

The detected spin is 1H; both the initial state and receiver are `L+` on 1H, and `decouple` is empty. Source settings are `offset=1000`, `sweep=2000`, `npoints=8192`, and `zerofill=65536`; the axis is ppm and inverted. Offset and sweep units are not given in this source. The FID is processed first with Gaussian apodisation parameter 5, then with the source's `bad-z1` parameter 10 and `bad-z2` parameter 40 operations (each call passes final argument 0). The shifted Fourier transform is plotted using its real part.

## Source limits

The source attributes the example to Damien Jeannerat and estimates minutes of calculation, but provides no DOI, explicit shim units, numerical spectrum, measured shim values, or relaxation model in this MATLAB file. The loaded MAT file supplies the spin-system numbers.
