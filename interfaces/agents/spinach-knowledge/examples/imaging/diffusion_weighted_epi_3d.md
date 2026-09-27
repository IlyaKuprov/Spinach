# examples/imaging/diffusion_weighted_epi_3d.m

- Signature: `diffusion_weighted_epi_3d()`

## Purpose

Simulates three-dimensional echo-planar imaging with diffusion using a brain phantom, then displays the acquired k-space data and a reconstructed image. Runtime is on the order of hours; the source notes that a Tesla V100 GPU can accelerate it, but the GPU setting is commented out.

## Physical / mathematical content

The example models a single ¹H spin at 5.9 T with zero chemical shift and diagonal T1/T2 relaxation, with both rates set to 1.0. It uses the `brain-medres` phantom for spatial R1, R2, and proton-density maps. The three-dimensional diffusion tensor field is isotropic: its diagonal entries are 2e-5 and its off-diagonal entries are zero. Slice selection, readout, phase encoding, and diffusion gradients are specified for the echo-planar sequence.

## Numerical / algorithmic content

The simulation uses a 201 × 201 image size, periodic spatial derivatives of order 3, and a 50-step Gaussian slice-selection pulse lasting 2.0e-4 s. It runs `imaging` with `epi_3d` and a 20 ms echo time. After acquisition, it halves the phase-encoding gradient amplitude for the plotting FOV and k-space extent, displays the fourth root of the signal for fringe visibility, applies squared-sine apodisation in both image dimensions, and reconstructs the real-space image with a shifted two-dimensional FFT.

## Implementation structure

The function creates a spin system and an unrestricted spherical-tensor Liouville-space basis, then sets pulse, gradient, phantom, relaxation, initial-state, detection, and diffusion parameters. It plots the phantom R1 map, runs the imaging simulation, and plots k-space beside the reconstructed image. Path tracing is disabled; `greedy` is enabled, while the optional `gpu` setting remains commented out.
