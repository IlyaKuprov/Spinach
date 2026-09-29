# examples/imaging/diffusion_weighted_epi_2d.m

Source: https://github.com/IlyaKuprov/Spinach/blob/main/examples/imaging/diffusion_weighted_epi_2d.m

## Experiment

Simulates a 2D Stejskal-Tanner spin-echo echo-planar diffusion-weighted acquisition, following Figure 1 of doi:10.1148/radiol.09090021 (https://doi.org/10.1148/radiol.09090021). The source estimates minutes of runtime and says a Tesla V100 can accelerate it; its GPU setting is commented out.

## Spin and image model

A single 1H spin at `sys.magnet=5.9` and zero scalar Zeeman shift is combined with the `brain-medres` phantom. The example takes slice 50 of the returned R1, R2, and proton-density volumes. R1/R2 superoperators are built with `rlx_t1_t2`; relaxation is set to the diagonal representation, zero equilibrium, and rates of 1.0 for both R1 and R2 (the source gives no units for those values). Initial state is `Lz`, receiver state is `L+`, and the coil profile is uniform. The spatial diffusion tensor is constant and isotropic: `dxx=dyy=1e-4`, with zero off-diagonal components. The image size is [101 105]; grid dimensions and point counts come from the selected phantom slice, with third-order periodic derivatives.

## Encoding and output

The `epi_2d` sequence uses readout and phase-encoding gradient amplitudes 5.3e-3 and 4.8e-3 T/m, both with configured duration 2e-3; diffusion-gradient amplitude is [1e-3 1e-3] T/m with duration 1e-2. The source does not annotate duration units. The displayed panels are the recorded image and the slice's R1 and R2 maps. After simulation, the example halves the phase-encoding amplitude for the plotting FOV calculation; this is a display adjustment, not a change to the simulated encoding. The cited DOI is retained from the source.
