# examples/imaging/dpfgse_signal_select.m

Source: https://github.com/IlyaKuprov/Spinach/blob/main/examples/imaging/dpfgse_signal_select.m

## Experiment and spin model

Demonstrates DPFGSE signal selection for a solution of GABA in water, with gradients and soft pulses modelled explicitly. The seven-spin system consists of seven 1H isotopes at `sys.magnet=5.9` Tesla. Chemical shifts are [3.00 3.00 1.88 1.88 2.28 2.28 4.80] ppm; scalar couplings are 7.36 Hz for spin pairs (1,3), (1,4), (2,3), and (2,4), and 7.58 Hz for (3,5), (3,6), (4,5), and (4,6). The basis uses `sphten-liouv`, `IK-2`, proximity level 1, and scalar-coupling connectivity. Path tracing and Krylov acceleration are disabled. The source estimates minutes of runtime and notes a Tesla V100 may speed it up, but GPU enablement is commented out.

## Selection and acquisition

The 1D domain is configured as 0.30 with 100 points and third-order periodic derivatives. Initial magnetisation and detection use `Lz` and `L+`, with uniform profiles; relaxation phantoms/operators are empty, and flow and diffusion are explicitly zero. Gradient amplitudes are [1e-3 1.5e-3], with `g_dur=1e-3`; the source does not annotate their units. Signal selection uses a ten-step Gaussian RF table with frequency 750, amplitude scale 2*pi*340, total configured duration 10e-3, and zero RF phase; RF units are not specified in the source. Maximum rank is 2.

Acquisition parameters are offset 800, sweep 1200, 512 points, and zero-fill to 2048; the plotted axis is configured in Hz and inverted. The source does not attach units directly to the offset/sweep values.

## Output

imaging runs `dpfgse_select`; the returned FID is exponentially apodised with parameter 6, Fourier-transformed after zero-filling and shifted, then plotted as the real spectrum. These are processing settings, not a reported experimental spectrum or measured signal.

Zero track elimination is explicitly enabled with `zte` in `sys.enable`.
