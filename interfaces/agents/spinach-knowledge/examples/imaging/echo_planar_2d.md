# examples/imaging/echo_planar_2d.m

## Purpose

Simulates 2D echo-planar imaging of a brain phantom with a single-proton model. It produces a recorded image and companion relaxation-map views.

## Spin and image model

The system is one `1H` spin at 5.9 T with zero chemical shift. Relaxation is `t1_t2`, with diagonal retention, zero equilibrium, and both configured rates set to 1.0; the basis uses `sphten-liouv` with no approximation. The source loads `brain-medres` relaxation and proton-density maps, takes slice 50 from each, and uses the first two entries of the phantom geometry and point-count arrays. The 2D image-size setting is `[101 105]`.

The proton-density map weights the initial `Lz` state; a uniform receive-coil map detects `L+`. The relaxation phantoms are paired with the `rlx_t1_t2` operators. Velocity fields are zero arrays and diffusion is set to zero.

## Gradient and acquisition settings

Readout and phase-encode gradient amplitudes are `5.3e-3` and `4.8e-3` T/m, respectively, as explicitly annotated in the source. Their duration fields are each set to `2e-3`, with no unit comment on those fields. The offset is zero and no decoupling is configured. Spatial differentiation uses `{'period',3}`.

## Output and caveats

The simulation calls `imaging` with `epi_2d`. Before plotting, the code halves `pe_grad_amp`, explaining that `G{1}` is effectively halved for the FOV calculation. The output figure shows the recorded image plus the `R1` and `R2` phantom maps. The source estimates seconds of runtime and says a Tesla V100 is faster, but the GPU-enable line is commented out; this is not a measured benchmark or an enabled-GPU run.

## Source

[MATLAB example](https://github.com/IlyaKuprov/Spinach/blob/main/examples/imaging/echo_planar_2d.m)
