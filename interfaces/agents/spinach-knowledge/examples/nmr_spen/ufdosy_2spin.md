# examples/nmr_spen/ufdosy_2spin.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_spen/ufdosy_2spin.m)

- Signature: `ufdosy_2spin()`
- Credits: Ludmilla Guduff, Jean-Nicolas Dumez, and Ilya Kuprov.

## Experiment and spin model

This is a simulated ultrafast diffusion-ordered NMR experiment for two coupled protons. It adds spatial flow and dipole-dipole (DD) and chemical-shift-anisotropy (CSA) relaxation to the diffusion encoding. The spin system is two 1H nuclei at 14.1 T, with chemical shifts 6.5 and 7.5 ppm and a 15 Hz scalar coupling. Each spin has a CSA tensor with eigenvalues [-10, -10, 20]; their listed Euler angles are [0, 0, 0] and [0, pi/2, 0]. The source lists spin coordinates [0, 0, 0] and [0, 0, 2.5] without stating coordinate units.

The calculation uses the full sphten-liouv basis, the NMR assumption, Redfield relaxation with a 1.0e-9 s correlation time, secular retention, and zero equilibrium. The relaxation operator and initial/detection states are supplied as spatial phantoms: relaxation uses a uniform spatial profile, the initial state is longitudinal 1H magnetisation, and detection is transverse 1H coherence. These are simulated states, not imported measurement data.

## Spatial encoding and acquisition

The modeled sample is 0.015 m long and represented by 3000 points; spatial derivatives use `parameters.deriv={'period',7}`. A uniform flow of 1.0e-4 m/s and diffusion coefficient 8.0e-10 m^2/s are set in the source.

The detected channel is 1H. The acquisition settings are 0.100 s acquisition time, 1.5e-6 s dwell, 128 points, 256 loops, a 0.52 T/m acquisition gradient, and a 3600 Hz offset. Spatial encoding uses 1000 pulse points, smfactor 0.1, Te=0.0015 s, Tau=0.0016 s, bandwidth 110000 Hz, gradient Ge=0.2535 T/m, and a smoothed chirp. The source calls `imaging` with `@spendosy`; it does not include measured data.

## Signal processing and scope

The returned signal is Fourier transformed with shifts along dimensions 1 and 2. The magnitude is displayed against a chemical-shift axis in ppm and a field-of-view axis in mm. The source describes the calculation as minutes on an NVIDIA Tesla A100 and much longer on CPU; this is a source estimate, not a reproduced timing. No external DOI, imported measurement, or experimental validation is specified by this example.
