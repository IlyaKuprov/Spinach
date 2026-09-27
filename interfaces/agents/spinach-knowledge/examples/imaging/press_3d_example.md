# examples/imaging/press_3d_example.m

- Signature: `press_3d_example()`

## Purpose

Calculates and plots a three-dimensional PRESS excitation profile with a tilted gradient system. The source notes that changing the pulse frequencies moves the hot spot through the sample; the estimated runtime is hours, faster with a Tesla V100 GPU.

## Model and sequence

The 3.0 T system contains two coupled protons with chemical shifts 1.5 and 3.7 and scalar coupling 10. It uses the `sphten-liouv` formalism with no basis approximation and enables the polyadic option. The sample is `[0.30 0.25 0.27] m`, sampled on `[108 90 111]` points; the image grid is `[101 103 105]`. The three slice-selection gradients are each `25 mT/m`, and the spatial-selection gradient is `5 mT/m` for `5e-4 s`.

Three phase-`pi/2` RF pulses use frequencies `{-30e3,+30e3,-50e3}`, amplitudes `2*pi*5000`, durations `{0.5e-4,1.0e-4,1.0e-4 s}`, and maximum ranks `{2,2,2}`. The gradient angles are `[pi/3 pi/4 pi/5]`.

## Output

The source runs `imaging(spin_system,@press_voxel_3d,parameters)` on a uniform-density phantom and plots the resulting excitation profile over the three field-of-view axes. It does not perform a spectral acquisition in this example.
