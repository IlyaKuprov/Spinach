# examples/imaging/diffusion_weighted_epi_3d.m

Source: https://github.com/IlyaKuprov/Spinach/blob/main/examples/imaging/diffusion_weighted_epi_3d.m

## Experiment

Simulates 3D echo-planar imaging with spatially uniform, isotropic diffusion and a brain phantom. The source estimates hours of runtime and says a Tesla V100 can accelerate it. `greedy` is enabled; the GPU setting is not—the code's GPU line is commented out.

## Spin, RF, and image model

The spin system contains one 1H spin at `sys.magnet=5.9` with zero scalar Zeeman shift, and uses diagonal T1/T2 relaxation with zero equilibrium. A `brain-medres` volume supplies R1, R2, and proton-density maps; the phantom dimensions and point counts configure the 3D grid. Initial and detection states are `Lz` and `L+`, with a uniform coil. The diffusion tensor is diagonal and uniform (`dxx=dyy=dzz=2e-5`); all cross terms are zero. The image size is [201 201].

Slice selection is represented explicitly by 50 Gaussian RF steps: total `pulse_time=2e-4`, `frequency table -5e3`, and `amplitude table 2*pi*7500*pulse_shape('gaussian',50)`. RF frequency/amplitude units and pulse-time units are not stated in the source. The slice-selection gradient is 32e-3 T/m for configured duration 1e-3; readout and phase-encoding gradients are 5.3e-3 and 4.8e-3 T/m, each with configured duration 4e-3. Diffusion-gradient amplitudes are [1e-3 1e-3 1e-3] T/m, echo-time parameter is 20e-3, and gradient angles are [pi/3 pi/4 pi/5]. Gradient-duration units are not annotated.

## Processing and caveats

The example plots the phantom R1 map, runs `epi_3d`, halves the phase-encoding amplitude for the subsequent FOV/k-space display, and shows the fourth root of the acquired signal to improve fringe visibility. It then applies squared-sine apodisation in both dimensions, computes a shifted 2D Fourier transform, takes its real part, and plots the reconstructed image. The reconstruction is thus the example's 2D display of data from the 3D acquisition, not a claim of a separately validated volumetric reconstruction. Runtime is source-estimated, not benchmarked here.
