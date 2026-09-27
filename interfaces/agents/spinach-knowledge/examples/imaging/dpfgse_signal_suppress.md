# examples/imaging/dpfgse_signal_suppress.m

- Signature: `dpfgse_signal_suppress()`

## Purpose

Simulates DPFGSE water suppression for a solution of GABA in water, with gradients and soft pulses implemented explicitly. The source estimates a simulation time of minutes, faster with a Tesla V100 GPU.

## Physical / mathematical content

- The seven-spin `1H` system has a 5.9 T magnet, chemical shifts of 3.00, 1.88, 2.28, and 4.80 ppm, and scalar couplings of 7.36 and 7.58 Hz.
- The sequence uses two gradient amplitudes, 1e-3 and 1.5e-3, with a duration of 1e-3 s. A ten-step Gaussian-shaped 180-degree pulse is specified for water.
- The sample uses a 0.30 spatial dimension with 100 grid points and a periodic third-derivative setting. Diffusion and flow are set to zero; no relaxation phantoms or operators are supplied.

## Numerical / algorithmic content

- The basis uses `sphten-liouv` formalism with `IK-2` approximation, proximity level 1, and scalar-coupling connectivity. Path tracing and Krylov acceleration are disabled.
- Acquisition uses a 800 Hz offset, 1200 Hz sweep, and 512 points. The FID is exponentially apodised with parameter 6, Fourier-transformed with 2048-point zero filling, shifted with `fftshift`, and plotted as the real spectrum in Hz with an inverted axis.

## Implementation structure

- Creates the spin system and basis, then sets initial `Lz` and detection `L+` states with uniform spatial phantoms.
- Specifies a ten-step pulse with frequency entries of 1220, Gaussian amplitude scaled by `2*pi*1700`, a total duration of 2e-3 s, zero phase, and maximum rank 2.
- Calls `imaging(spin_system,@dpfgse_suppress,parameters)`, processes the resulting FID, and plots the spectrum. GPU enablement appears only as a commented-out setting.
