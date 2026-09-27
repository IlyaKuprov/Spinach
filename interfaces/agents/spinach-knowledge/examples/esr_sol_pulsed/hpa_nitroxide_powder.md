# examples/esr_sol_pulsed/hpa_nitroxide_powder.m

- Signature: `hpa_nitroxide_powder()`

## Purpose

Simulates the powder-averaged pulse-acquire W-band ESR spectrum of an electron–¹⁴N nitroxide radical at 3.5 T. The script models acquisition and Fourier processing; it does not define an excitation pulse sequence. Calculation time: seconds.

## Physical model

The two-spin system contains an electron (E) and ¹⁴N. The electron g tensor and electron–nitrogen hyperfine tensor are anisotropic and include off-diagonal components. The spin system uses a secular diagonal relaxation model with a 5×10⁷ s⁻¹ damping rate and zero equilibrium state.

## Simulation and processing

The calculation uses the full-sphten Liouville-space basis without a basis approximation and disables trajectory-level SSR. The initial state and detection coil are both the electron `L+` operator. The acquisition uses a 1 GHz sweep, 128 points, an offset of −2×10⁸ (in the script's frequency units), and the `rep_2ang_6400pts_sph` powder grid; the spectrum axis is labelled in GHz in the lab frame and inverted. After powder averaging, the FID is apodised with `crisp`, zero-filled to 512 points, Fourier transformed, and the real spectrum is plotted.
