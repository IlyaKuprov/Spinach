# examples/nmr_zerofield/zero_field_benzene.m

## Experiment represented

This is a simulated zero-field NMR spectrum for benzene with one 13C nucleus. The source says the example is set to reproduce Figure 2 of [the cited Journal of Magnetic Resonance paper](https://doi.org/10.1016/j.jmr.2017.08.016). The spectrum is calculated from the source-defined spin system rather than imported from measurement.

## Spin model

The ordered spin list is six 1H nuclei followed by one 13C nucleus, and sys.magnet is zero. The source explicitly assigns scalar-coupling values (Hz): C7 couples to H1–H6 by 158.363, 1.136, 7.609, −1.285, 7.609, and 1.136, respectively. The proton-pair couplings are H1–H2 7.534, H1–H3 1.381, H1–H4 0.658, H1–H5 1.381, H1–H6 7.534; H2–H3 7.543, H2–H4 1.382, H2–H5 0.660, H2–H6 1.384; H3–H4 7.543, H3–H5 1.387, H3–H6 0.660; H4–H5 7.543, H4–H6 1.382; and H5–H6 7.543. The basis is zeeman-hilb with approximation none. No temperature, gradient, or chirp is specified.

## Acquisition and processing

The simulation uses a 400 Hz sweep, 8162 points, and 16586 zero-fill points; offset is zero, 1H is the selected channel, axis units are Hz, and inversion is disabled. The nominal flip angle is π/2 with uniaxial detection. The source builds the system and basis, then calls liquid with zerofield in the lab frame. It subtracts the FID mean, applies exponential apodisation with parameter 6, computes the shifted FFT at the zero-fill length, and plots the real spectrum. The header calls this a seconds-scale calculation; that is not a measured runtime here.

## Source

[MATLAB example](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_zerofield/zero_field_benzene.m)