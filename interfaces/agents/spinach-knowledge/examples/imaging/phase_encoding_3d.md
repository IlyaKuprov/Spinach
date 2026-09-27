# examples/imaging/phase_encoding_3d.m

- Signature: `phase_encoding_3d()`

## Purpose

Simulates 3D slice selection followed by phase-encoded imaging of the selected slice. The example uses the library brain-medres phantom and reconstructs an image from the acquired signal. The source estimates minutes of runtime, with a Tesla V100 GPU reported as faster.

## Model and sequence

The spin system is one `1H` at 5.9 T, with zero chemical shift, diagonal T1/T2 relaxation, rates `r1 = r2 = 1`, and zero equilibrium. The example obtains the R1, R2, and proton-density maps and dimensions from `phantoms('brain-medres')`; diffusion and all three flow components are set to zero.

A 50-step Gaussian RF pulse (total duration `2.0e-4 s`, frequency `-5 kHz`, peak scale `2*pi*7500`, phase `pi/2`) performs slice selection. The imaging grid is `[129 129]`; slice-selection, phase-encoding, and readout gradients are each `32 mT/m`, with durations `1.0e-4`, `2.0e-4`, and `3.0e-4 s`, respectively. The echo time is `20 ms`, and the gradient angles are `[pi/3 pi/4 pi/5]`.

## Computation and output

The sequence is run with `imaging(spin_system,@phase_enc_3d,parameters)`. The example displays the k-space data, applies square-sine apodisation in both dimensions, computes a shifted 2D Fourier transform, and displays the real-space image. Path tracing and Krylov propagation are disabled in the source; the greedy option is enabled.
