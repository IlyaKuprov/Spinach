# examples/nmr_metabol/molecule_a.m

- Signature: `molecule_a()`
- Source: [examples/nmr_metabol/molecule_a.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_metabol/molecule_a.m)

## What it models and loads

This example simulates a one-dimensional liquid-state `1H` NMR spectrum for a molecule identified in the source as a GISSMO-database entry. The source comment estimates calculation time as seconds; that is not a timing measurement from this review. The wrapper imports `molecule_a.xml` with `gissmo2spinach('molecule_a.xml',1)` to obtain `sys` and `inter`. That XML is an input to the spin-system simulation, not an acquired FID or spectrum loaded for plotting. This wrapper does not say whether the XML parameters were experimentally measured, calculated, or curated.

## Spin system and acquisition setup

- The basis uses `sphten-liouv` formalism, `IK-2` approximation, `scalar_couplings` connectivity, and proximity level `1`.
- The observed spins are `{'1H'}`; both the initial state and detection coil are set to `state(spin_system,'L+','1H')`. `parameters.decouple={}` specifies no decoupling entries.
- The source sets offset `3500`, sweep `3000`, `4096` acquisition points, and `16536` zero-filled points. It labels the plotted axis `ppm` and sets `invert_axis=1`. The wrapper does not state units for offset, sweep, or the Gaussian parameter below; these numbers are the literal settings in the code, not inferred Hz or ppm values.
- No relaxation parameters are assigned in this wrapper.

Acquisition is delegated through `liquid(spin_system,@acquire,parameters,'nmr')`. The wrapper does not show the helper's internal pulse, phase, gradient, or receiver program, so those details are not asserted here.

## Processing and output

The returned FID is apodised with `{'gauss',10}`, Fourier-transformed as `fftshift(fft(fid,parameters.zerofill))`, and displayed as `real(spectrum)` with `plot_1d`. The zero-argument function does not declare a returned value or write an output spectrum file in this source. No experimental spectrum is loaded or compared in the wrapper.
