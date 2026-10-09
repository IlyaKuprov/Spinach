# examples/esr_sol_pulsed/hpa_gd_dota_powder.m

- MATLAB implementation: [examples/esr_sol_pulsed/hpa_gd_dota_powder.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/esr_sol_pulsed/hpa_gd_dota_powder.m)

[MATLAB example](../../../../../examples/esr_sol_pulsed/hpa_gd_dota_powder.m) · [acquire sequence helper](../../../../../experiments/acquire.m)

- Signature: `hpa_gd_dota_powder()`

## Aim and spin system

This example calculates a powder-averaged pulsed ESR spectrum for a Gd(III)–DOTA complex; its source describes an ideal pulse, W-band operation, and a numerical second-order rotating-frame transformation. The model uses one `E8` electron spin at 9.40 T, isotropic Zeeman scalar `1.9918`, and an axial ZFS tensor with eigenvalues `[0.57e9, 0.57e9, -2*0.57e9]/3` and Euler angles `[0 0 0]`. The Zeeman Liouville basis is exact (`approximation='none'`), and trajectory-level SSR algorithms are disabled.

## Acquisition protocol

The script sets both initial state and detection operator to `L+` on `E8`, the source's representation of the ideal-pulse transverse coherence and detected channel; it supplies these to `acquire` rather than simulating an RF pulse explicitly. No spins are decoupled. The fixed acquisition settings are offset `1.5e9` (unit not annotated), sweep width `6e9` Hz, 4096 points, and zero-fill to 16384. It uses the `rep_2ang_12800pts_sph` spherical powder grid, second-order rotating frames for `E8`, no derivative, and an inverted axis labelled `GHz-labframe`. The `acquire` helper defines sweep width in Hz and evolves the specified initial state to form the FID.

## Observable and output

The powder-averaged FID is exponentially apodised with parameter 6, Fourier transformed, and its real spectrum is plotted. No output file is written by the example. The source estimates a calculation time of minutes; its stated ideal-pulse and second-order rotating-frame assumptions delimit the simulation.
