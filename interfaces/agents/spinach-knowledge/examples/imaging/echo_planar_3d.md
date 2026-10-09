# examples/imaging/echo_planar_3d.m

## Purpose

Demonstrates 3D slice selection followed by echo-planar acquisition on the `brain-medres` phantom. Despite the 3D sequence and volumetric phantom, the example displays k-space for a slice and reconstructs a 2D image using a 2D FFT; it is not a full 3D image reconstruction.

## Spin and image model

The model is one `1H` spin at 5.9 T with zero chemical shift, `t1_t2` relaxation, diagonal retention, zero equilibrium, and both rate settings equal to 1. The basis is `sphten-liouv` with no approximation. The phantom call supplies full 3D `R1`, `R2`, and proton-density maps, dimensions, and point counts. The proton-density phantom weights the initial `Lz` state; a uniform coil vector detects `L+`, and `rlx_t1_t2` supplies the relaxation operators. Flow fields and diffusion are set to zero. The `R1` volume is plotted before the sequence.

## Slice-select RF and gradients

The slice-select pulse is a 50-step Gaussian pulse with phase `pi/2`, frequency entries `-5e3`, amplitude scaled by `2*pi*7500`, and total duration parameter `2.0e-4` divided equally among the steps. These RF-table values and durations have no unit comments in the source. The image-size array is `[201 201]`. Slice-select, readout, and phase-encode gradient amplitudes are `32.0e-3`, `5.3e-3`, and `4.8e-3` T/m, as labelled in the source. Readout and phase-encode duration fields are both `4e-3`; the echo-time field is `20e-3`. The gradient-angle array is `[pi/3 pi/4 pi/5]` and has no separate unit annotation. Spatial differentiation uses `{'period',3}`. The source actively sets `sys.enable={'zte','greedy'}`; the nearby `gpu` text is a comment, not an enabled GPU option.

## Output and caveats

The sequence call is `imaging(...,@epi_3d,...)`. The code halves `pe_grad_amp` for the plotted FOV and k-space extent because `G{1}` is effectively halved. It displays `fid.^(1/4)` as the k-space view (the source says this improves fringe visibility), then applies square-sine apodisation in both dimensions and forms a real-space slice with `fft2`. The source estimates hours of runtime and says a Tesla V100 is faster; this file does not enable GPU execution.

## Source

[MATLAB example](https://github.com/IlyaKuprov/Spinach/blob/main/examples/imaging/echo_planar_3d.m)
