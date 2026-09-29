# examples/nmr_zerofield/field_drop_acetonitrile.m

## Experiment represented

This is a simulated zero-field NMR field-drop experiment for acetonitrile. The source comment describes computing the exact thermal equilibrium state and propagating it through a time-dependent field drop; it identifies the plotted target as Figure 7 of [the cited Journal of Magnetic Resonance paper](https://doi.org/10.1016/j.jmr.2017.08.016). This is a Spinach calculation, not a file of imported measured data.

## Spin model and preparation

The ordered spin list is three 1H nuclei, two 13C nuclei, and one 15N nucleus. The polariser setting is 2.0 T, and the source sets the temperature to 298 K. Scalar-coupling entries (Hz) are: each of H1–H3 to C4, 136.200; each of H1–H3 to C5, −9.924; each of H1–H3 to N6, −1.688; C4–C5, 57.010; C4–N6, 2.822; and C5–N6, −17.419. The basis uses the full sphten-liouv formalism with approximation none and S3 permutation symmetry on spins 1–3.

## Field-drop sequence and computed spectrum

The calculation sets a 700 Hz sweep, 4096 acquired points, 16586 zero-fill points, zero offset, and 1H observation in Hz without axis inversion. Its drop settings are 0.1 mT final field, 10 Hz drop rate, 0.5 s drop time, and 100 drop points; the nominal flip angle is π/2 and detection is uniaxial. The script passes this configuration to liquid with zulf_abrupt in the lab frame. It subtracts the FID mean, applies exponential apodisation with parameter 6, computes fftshift(fft(fid, zerofill)), and plots the real spectrum. The script does not specify a gradient or chirp schedule.

The source header describes the calculation time as “seconds”; that is source commentary, not a timing measurement made here. No run or validation was performed for this knowledge-page update.

## Source

[MATLAB example](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_zerofield/field_drop_acetonitrile.m)