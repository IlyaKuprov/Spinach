# examples/imaging/diffusion_weighted_epi_2d.m

- Signature: `diffusion_weighted_epi_2d()`

## Purpose

2D echo planar imaging example in the presence of isotropic diffusion. Stejskal-Tanner SE echo planar diffusion-weighted pulse sequence from Figure 1 in (https://doi.org/10.1148/radiol.09090021). Simulation time: minutes, faster with a Tesla V100 GPU.

## Physical / mathematical content

Simulates a 2D Stejskal–Tanner spin-echo EPI acquisition with isotropic diffusion. The spatial diffusion tensor has dxx=dyy=1e-4 and dxy=dyx=0 throughout the grid. A brain phantom supplies proton density and spatially varying R1 and R2 relaxation; the initial spin state is 1H Lz and detection uses 1H L+ with a uniform coil profile.

## Numerical / algorithmic content

Creates a single-1H spin system at 5.9 T with zero chemical shift, diagonal T1/T2 relaxation, zero equilibrium, and R1 and R2 rates of 1.0. Uses the sphten-liouv formalism without basis approximation, disables path tracing, and configures periodic third-order spatial derivatives. The image size is [101 105]; readout and phase-encoding gradients are 5.3e-3 and 4.8e-3 T/m, each lasting 2e-3 s, and the diffusion-gradient amplitudes are [1e-3 1e-3] T/m with duration 1e-2 s. The simulation calls imaging with epi_2d, then halves the phase-encoding gradient amplitude for field-of-view calculation during plotting.

## Implementation structure

Builds the spin system and basis, obtains T1/T2 relaxation superoperators, and takes slice 50 of the brain-medres R1, R2, and proton-density phantoms. It configures the phantom geometry, relaxation maps, spin states, coil, and diffusion field before running the EPI simulation. Three panels display the recorded image and the R1 and R2 phantoms. GPU enablement is present only as a commented-out setting.

- Figure 1 in (https://doi.org/10.1148/radiol.09090021).
