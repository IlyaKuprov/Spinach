# examples/parahydrogen/altadena_propanal.m

- Signature: altadena_propanal()
- Source: [examples/parahydrogen/altadena_propanal.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/parahydrogen/altadena_propanal.m)

## Model and assumptions

This is a Spinach simulation of ALTADENA polarisation following parahydrogenation of acrolein to propanal. Para-hydrogen carries correlated proton-pair spin order into the product; under the idealised ALTADENA picture, adiabatic transport from the low-field reaction region to the high-field detection region converts that order into observable product polarisation. The script represents the product spin system and its chosen initial density operator; it does not simulate the chemical addition, a field ramp, or a measured enhancement. Its header explicitly assumes perfectly adiabatic transfer and omits isotropic mixing at low field, so the plotted signal is conditional on those simplifications.

The model has six spin-1/2 protons at 7.05 T and uses the `sphten-liouv` formalism without basis approximation.  The three equivalent sites assigned 1.11 ppm and the two sites at 2.46 ppm are grouped with S3 and S2 symmetry, respectively; the remaining site is assigned 9.79 ppm. Scalar couplings are 7.3 Hz from each of the first three spins to each of spins 4 and 5, and 1.4 Hz from spins 4 and 5 to spin 6. The source expresses the chemical shifts as values in `inter.zeeman.scalar`; the acquisition axis is configured in ppm.

## Initial order and computed signal

The initial operator weights the two-spin longitudinal term for spins 1 and 4 by 1.0, the spin-1 longitudinal term by −0.5, and the spin-4 term by +0.5. A small 1H pulse of `pi/100` (1.8 degrees) converts part of the resulting order into a detectable transverse signal. `liquid` with `hp_acquire` computes the simulated free-induction signal; exponential apodisation (parameter 6), Fourier transformation, zero filling to 8192 points, and plotting produce the spectrum. Acquisition is set to a 500 Hz offset, 1000 Hz sweep, and 1024 acquired points, with the display axis in ppm and reversed orientation. The source comments estimate seconds of calculation time; this is a code comment, not a benchmark performed here.

Source authors in the existing entry are Ronghui Zhou (hui@ufl.edu) and Ilya Kuprov (ilya.kuprov@weizmann.ac.il). No DOI or published hyperpolarisation measurement is given in this example. The result should therefore be read as an idealised simulated ALTADENA spectrum, not experimental evidence for a measured polarisation level.
