# examples/nmr_zerofield/small_field_formic_acid.m

## Experiment represented and source-label discrepancy

The filename and function name identify this as the small-field formic-acid example, while its header comment instead says “15N pyridine” and says it is set to reproduce Figure 2 of [the cited Physical Review Letters paper](https://doi.org/10.1103/PhysRevLett.107.107601). The executable spin model is unambiguous about its channels: it contains one 1H and one 13C, not 15N. The script computes a simulated signal; it does not import measured data.

## Spin model and acquisition

The code sets sys.magnet to 1.76e-7 T and defines one 1H–13C scalar coupling of 221 Hz. It does not set a temperature. The basis is zeeman-hilb with approximation none. Acquisition uses a 600 Hz sweep, 8196 points, and 16384 zero-fill points; offset is zero, the selected detection channel is 1H, the axis is in Hz without inversion, and the nominal flip angle is π/2 with uniaxial detection. No gradient or chirp schedule is specified.

Spinach builds the spin system and basis, then liquid propagates it with zerofield in the lab frame. The FID is mean-subtracted and exponentially apodised with parameter 12, followed by a shifted FFT using the zero-fill length and a plot of the real spectrum. The source header estimates a calculation time of seconds.

## Source

[MATLAB example](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_zerofield/small_field_formic_acid.m)