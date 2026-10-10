# examples/nmr_spen/psycosy_salsalate.m

- Function: `psycosy_salsalate()`
- Source: [examples/nmr_spen/psycosy_salsalate.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_spen/psycosy_salsalate.m)

This is a simulated PSYCOSY experiment for a five-proton model of one salsalate ring. The source supplies a 14.1 T magnet, five 1H shifts (8.14, 7.44, 7.71, 7.28, and 7.60 ppm), and scalar couplings of 7.9, 7.5, 8.1, 1.6, and 1.2 Hz; it also assigns zero couplings for the listed (5,2) and (5,5) entries. These are simulation inputs, not imported experimental measurements.

The spatial sample is 15 mm long and represented by 100 points, with the period-3 spatial derivative setting. Flow is zero and the diffusion coefficient is zero, so this calculation does not model diffusion attenuation. The initial and receive phantoms are uniform; the spin states are 1H longitudinal magnetisation and 1H transverse detection. Relaxation is represented by the relaxation operator with a zero spatial relaxation phantom.

The script calls Spinach imaging with the PSYCOSY sequence. It sets a 1H, 720 Hz sweep, offset 4620, and a 512 by 512 acquisition, zero-filled to 1024 by 1024. The sequence inputs include 110 ms mixing and a 0.01 T/m gradient. Its saltire-chirp inputs are 20 degrees, 15 ms pulse duration, 50 ms gradient duration, 10 kHz sweep, 250 points, and smoothing factor 20. The source passes these parameters to the sequence; it does not spell out the internal pulse/gradient waveform or a detailed pathway derivation, so none is inferred here.

Each FID dimension receives square-sine apodisation; a two-dimensional FFT is shifted and plotted as a positive magnitude spectrum. The model is configured with the sphten-liouv formalism, IK-2 approximation, proximity level 1, scalar-coupling connectivity, greedy algorithm, proximity cutoff 4.0, and merge dimension 500. These settings describe this example only; the file contains no comparison with measured data or stated accuracy bound. The source comment's machine-time estimate is not repeated as a validated runtime.

Zero track elimination is explicitly enabled with `zte` in `sys.enable`.
