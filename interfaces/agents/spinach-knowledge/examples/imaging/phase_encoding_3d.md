# examples/imaging/phase_encoding_3d.m

Three-dimensional slice selection followed by phase-encoded imaging of the selected slice. The source estimates minutes of simulation time and says a Tesla V100 is faster; it credits Ahmed Allami and Ilya Kuprov. [Source](https://github.com/IlyaKuprov/Spinach/blob/main/examples/imaging/phase_encoding_3d.m)

## Model and sequence

The system is a single `1H` spin with `sys.magnet=5.9` and zero scalar chemical shift (the example does not annotate the field-value unit). It uses `t1_t2` relaxation with `R1=1.0`, `R2=1.0`, zero equilibrium, and diagonal retention. The `sphten-liouv` basis has no approximation. The source disables path tracing and Krylov propagation; it enables `zte` and `greedy` while the adjacent `'gpu'` text is a comment, not the active setting. The source comment says a GPU is needed.

`phantoms('brain-medres')` supplies the `R1`, `R2` and proton-density maps plus sample dimensions and point counts. Relaxation operators come from `rlx_t1_t2`; the `R1`/`R2` maps are used as relaxation phantoms and proton density as the initial spatial phantom. The initial state is `Lz`, the uniform coil uses `L+`, and all three flow fields are zero with diffusion `0`. Grid differentiation is `{'period',3}`.

Slice selection is set by a phase of `pi/2` and 50 pulse steps over `2.0e-4` total duration. The source builds a 50-step frequency list at `-5e3` and a Gaussian-shaped amplitude list from `2*pi*7500*pulse_shape('gaussian',50)`; it does not annotate units for these RF values or pulse duration. Image size is `[129 129]`. Slice-select, phase-encode and readout gradient amplitudes are each `32.0e-3 T/m`; their duration values are `1.0e-4`, `2.0e-4` and `3.0e-4`, with no units stated. The echo value is `20e-3` and gradient angles are `[pi/3 pi/4 pi/5]`.

## Output and interpretation

The script first plots the 3D `R1` phantom, then calls `imaging(spin_system,@phase_enc_3d,parameters)` for slice data. It displays the raw data as `fid.^(1/4)` in k-space, applies square-sine apodisation in both dimensions, computes `real(fftshift(fft2(ifftshift(fid))))`, and displays the reconstructed 2D slice. The transformed display is a 2D slice result, not a 3D Fourier reconstruction.
