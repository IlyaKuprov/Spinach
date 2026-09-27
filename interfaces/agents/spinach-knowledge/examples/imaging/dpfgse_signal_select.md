# examples/imaging/dpfgse_signal_select.m

- Signature: `dpfgse_signal_select()`

## Purpose

Simulate DPFGSE signal selection for a solution of GABA in water, with gradients and soft pulses represented explicitly. The source estimates a runtime of minutes, faster with a Tesla V100 GPU.

## Physical / mathematical content

- The seven-spin `1H` system uses a 5.9 T magnet, chemical shifts of 3.00, 3.00, 1.88, 1.88, 2.28, 2.28, and 4.80 ppm, and scalar couplings of 7.36 and 7.58 Hz between the specified spin pairs.
- The simulation uses a 0.30-unit spatial domain with 100 points and `parameters.deriv={'period',3}`. Flow and diffusion are set to zero; the relaxation phantom and operator arrays (`parameters.rlx_ph` and `parameters.rlx_op`) are empty.
- Uniform initial-state and coil phantoms use `Lz` and `L+` proton states, respectively.

## Numerical / algorithmic content

- The basis uses `sphten-liouv` formalism, `IK-2` approximation, proximity level 1, and scalar-coupling connectivity. Path tracing and Krylov acceleration are disabled.
- Signal selection uses gradient amplitudes of 1e-3 and 1.5e-3 with a duration of 1e-3, and a ten-step Gaussian RF pulse with frequency 750, amplitude scale `2*pi*340`, and total duration 10e-3. The maximum rank is 2 and the RF phase is zero.
- Acquisition specifies a proton offset of 800, sweep of 1200, and 512 points. The resulting FID is exponentially apodised with parameter 6, Fourier-transformed with zero filling to 2048 points, shifted with `fftshift`, and plotted as a real spectrum on an inverted Hz axis.

## Implementation structure

- Construct the spin system and basis, then pass the sequence, spatial-grid, and state parameters to `imaging(spin_system,@dpfgse_select,parameters)`.
- Process the returned FID with `apodisation` and `fftshift(fft(...))`, then plot it with `plot_1d`.
- GPU enabling appears only as a commented-out setting; the function does not enable it.
