# examples/esr_sol_pulsed/hpa_nitroxide_powder.m

- MATLAB implementation: [examples/esr_sol_pulsed/hpa_nitroxide_powder.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/esr_sol_pulsed/hpa_nitroxide_powder.m)

[MATLAB example](../../../../../examples/esr_sol_pulsed/hpa_nitroxide_powder.m) · [acquire sequence helper](../../../../../experiments/acquire.m)

- Signature: `hpa_nitroxide_powder()`

## Aim and spin system

The source describes a powder-averaged pulse-acquire W-band Fourier ESR spectrum for a nitroxide radical and assumes an ideal pulse. The model has an electron (`E`) and `14N` at 3.5 T. The electron g tensor is `diag(2.01045, 2.00641, 2.00211)`; the electron–nitrogen coupling tensor is supplied as `1e7 * [1.2356 0 0.6322; 0 1.1266 0; 0.6322 0 8.2230]` (the source does not annotate the matrix unit). The source also specifies damping relaxation (`inter.relaxation={'damp'}`), retains diagonal relaxation terms, sets equilibrium to zero, and sets `damp_rate=5e7` (the source does not annotate this field's unit). It uses an exact spherical-tensor Liouville basis and disables trajectory-level SSR algorithms.

## Acquisition protocol

Both the initial state and coil operator are `L+` on the electron. The code passes these to `acquire`, so the ideal pulse is represented by the chosen transverse initial coherence; no explicit pulse waveform, duration, or pulse-power scan is simulated. No spins are decoupled. Fixed acquisition settings are offset `-2e8`, sweep width `1e9` Hz, 128 points, zero-fill to 512, and the `rep_2ang_6400pts_sph` powder grid. The derivative is disabled, the axis is inverted, and the display is labelled `GHz-labframe`. The `acquire` helper documents sweep width in Hz.

## Observable and output

The powder-averaged FID is apodised with the `crisp` window, Fourier transformed, and the real spectrum is plotted. The example does not save a data or figure file and estimates a run time of seconds. This is a simulated ideal-pulse acquisition with the specified damping model, not a measured spectrum.
