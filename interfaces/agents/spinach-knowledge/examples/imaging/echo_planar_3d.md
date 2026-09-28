# examples/imaging/echo_planar_3d.m

- Signature: `echo_planar_3d()`

## Purpose

Simulate 3D slice selection followed by a three-dimensional echo-planar imaging sequence using a brain phantom, then display a reconstructed 2D slice. The source notes a simulation time of hours, faster with a Tesla V100 GPU; GPU execution is not enabled in this file.

## Physical / mathematical content

- Models a single `1H` spin at 5.9 T with zero chemical shift and diagonal T1/T2 relaxation, using rates of 1 for both R1 and R2. The `brain-medres` phantom supplies spatial R1, R2, and proton-density maps.
- Uses a 50-step Gaussian slice-selection pulse of duration `2.0e-4` s, with RF frequency `-5e3`, amplitude scale `2*pi*7500`, and phase `pi/2`. Slice-selection, readout, and phase-encoding gradient amplitudes are `32.0e-3`, `5.3e-3`, and `4.8e-3` T/m; gradient angles are `[pi/3 pi/4 pi/5]`.
- Sets diffusion and flow to zero. The initial spin state is `Lz`, and detection uses `L+` with a uniform coil phantom.

## Numerical / algorithmic content

- Creates a spin system in the `sphten-liouv` formalism with no basis approximation, then runs `imaging(spin_system,@epi_3d,parameters)` with image size `[201 201]`, 4 ms readout and phase-encoding gradient durations, a 20 ms echo time, and periodic third-order spatial derivatives.
- Displays the 3D R1 map, plots the acquired slice's k-space data after a fourth-root display scaling, applies `sqsin` apodisation in both dimensions, and reconstructs the real-valued image with a shifted 2D FFT. The phase-encoding gradient amplitude is halved after simulation for the plotting field of view and k-space extent.

## Implementation structure

- Defines the spin system, relaxation model, basis, RF pulse, and sequence parameters.
- Loads the `brain-medres` phantom and supplies relaxation, initial-state, and detection phantoms to `epi_3d` through `imaging`.
- Plots the phantom and the simulated slice in k-space and real space.
