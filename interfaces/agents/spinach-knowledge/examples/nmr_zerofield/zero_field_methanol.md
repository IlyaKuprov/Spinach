# examples/nmr_zerofield/zero_field_methanol.m

- Signature: `zero_field_methanol()`
- Source: [`examples/nmr_zerofield/zero_field_methanol.m`](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_zerofield/zero_field_methanol.m)

## Purpose

This is a simulated zero-field NMR example for 13C-labelled methanol, set up to reproduce Figure 1 of [Chemical Physics Letters, DOI 10.1016/j.cplett.2013.06.042](https://doi.org/10.1016/j.cplett.2013.06.042). The reference is the example's stated target, not a claim that this copy has been run or independently matched to the figure.

## Spin and field model

The four-spin system contains three 1H spins and one 13C spin. The field is set to zero, and the scalar 1H–13C couplings from each proton to carbon are each 141.0 Hz. The source selects the spherical-tensor Liouville formalism with no approximation and uses S3 permutation symmetry for the three proton spins. No spatial coordinates, gradient, chirp, field-drop schedule, or external measured trace is specified.

## Acquisition and processing

The zero-field signal is calculated with `liquid(...,@zerofield,...,'labframe')`. The acquisition settings specify a 700 Hz sweep, 4096 points, 16384-point zero filling, a 1H channel, a π/2 flip angle, uniaxial detection, and a zero offset; the plotted axis is in Hz. The code subtracts the FID mean, applies exponential apodisation with parameter 6, Fourier transforms and shifts the result, then plots its real part. The source comment estimates calculation time as seconds; that is a source note, not a runtime measured here.

This is a zero-field NMR simulation, not an SPEN or ultrafast-DOSY experiment: there is no encoded spatial dimension or gradient schedule, and the source does not define a multiple-quantum selection step. It imports no experimental data.
