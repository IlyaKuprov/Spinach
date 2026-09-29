# examples/nmr_spen/ufdosycosy_2spin.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_spen/ufdosycosy_2spin.m)

- Signature: `ufdosycosy_2spin()`
- Credits: Ludmilla Guduff and Jean-Nicolas Dumez.

## Experiment and spin model

This example simulates an ultrafast three-dimensional DOSY-COSY data set for a coupled two-proton system. Both isotopes are 1H at 14.1 T; the chemical shifts are 2.0 and -1.3 ppm and the scalar coupling is 15 Hz. It uses the full sphten-liouv basis and an NMR assumption. The source has empty relaxation phantom and operator lists, so it does not set up a relaxation contribution in this simulation. Its uniform phantom starts from longitudinal 1H magnetisation and detects transverse 1H coherence. No experimental data are loaded.

## Spatial encoding and sequence parameters

The sample length is 0.015 m with 3000 spatial points and derivative setting `parameters.deriv={'period',7}`. Flow is set to zero at all points; the diffusion coefficient is 4.0e-10 m^2/s.

The detected channel is 1H. Acquisition settings are 0.100 s acquisition time, 1.5e-6 s dwell, sweep parameter 1302, 128 points in dimension 1, 64 in dimension 2, 256 loops, acquisition gradient 0.52 T/m, and offset 600 Hz. Under the third-dimension comment, the source assigns the same sweep value, 1302 Hz. Encoding uses 1000 pulse points, smfactor 0.1, Te=0.0015 s, Tau=0.0016 s, bandwidth 110000 Hz, Ge=0.2535 T/m, and a smoothed chirp.

The source also sets diffusion scale dscale=1.0e-10, FOV limits -0.0041 and +0.0041 (units are not stated there), Hamming apodisation, a sine window, and the keeler_corr fitting model. These are simulation/model settings, not evidence of a fit to imported measurements. The source calls `imaging` with `@spendosycosy`.

## Reconstruction and source-defined scope

The returned array is Fourier transformed, with shifts, along all three dimensions. The plotted quantity is the magnitude raised to the one-half power, displayed as a volume over limits [-1, 1] for each plotted coordinate. The source does not name physical units for those three plotted limits, so they are not interpreted here. It estimates hours on an NVIDIA Tesla A100 and much longer on CPU; this is not a measured runtime. The file specifies a three-dimensional processing path but does not explicitly label the physical meaning of each transformed axis beyond its DOSY-COSY context and the listed acquisition parameters.
