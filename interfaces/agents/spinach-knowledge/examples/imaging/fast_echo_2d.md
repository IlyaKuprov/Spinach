# examples/imaging/fast_echo_2d.m

## Purpose

Runs the 2D fast spin-echo brain-imaging example using a single slice of the `brain-medres` phantom. “Fast” describes the experiment-duration intent in the source comment; the same header estimates simulation time in hours.

## Spin and image model

The model is one `1H` spin at 5.9 T with zero chemical shift. It uses `t1_t2` relaxation, diagonal retention, zero equilibrium, and rate settings of 1 for both `R1` and `R2`; the basis is `sphten-liouv` with no approximation. The code selects slice 50 from the `R1`, `R2`, and proton-density maps, uses the first two geometry and point-count entries, and sets image size to `[101 105]`.

The proton-density map weights the initial `Lz` state, and a uniform receive-coil phantom detects `L+`. The relaxation maps are paired with the `rlx_t1_t2` operators. The source does not configure flow or diffusion parameters in this example.

## Gradient and acquisition settings

Readout and phase-encode amplitudes are `5.3e-3` and `4.8e-3` T/m, explicitly labelled in the source. Their duration fields are `2e-3` and `1e-3`, respectively; the source does not annotate duration units. The offset is zero, decoupling is empty, and spatial differentiation uses `{'period',3}`.

## Output and caveats

The example calls `imaging` with `fse`, then plots the recorded image beside the selected `R1` and `R2` phantom maps. The source estimates hours of runtime and says a Tesla V100 is faster; its GPU-enable line is commented out, so this file does not enable GPU execution. The timing is an estimate in the source, not a benchmark.

## Source

[MATLAB example](https://github.com/IlyaKuprov/Spinach/blob/main/examples/imaging/fast_echo_2d.m)
