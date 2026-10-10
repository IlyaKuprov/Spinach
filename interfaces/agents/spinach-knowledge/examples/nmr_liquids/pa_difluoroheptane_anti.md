# examples/nmr_liquids/pa_difluoroheptane_anti.m

- Signature: `pa_difluoroheptane_anti()`

## Purpose

Simulates a one-dimensional pulse-acquire 1H NMR spectrum for anti-3,5-difluoroheptane. The manually specified 23-spin isotope list contains seven 12C, fourteen 1H, and two 19F spins. The source cites [DOI 10.1021/acs.joc.4c00670](https://doi.org/10.1021/acs.joc.4c00670) for the manual basis construction; its header estimates minutes for calculation time and says it is faster with a GPU.

## Spin model and acquisition

The field setting is `11.7464`. The source assigns proton chemical-shift values `1.0092` and `4.6834`. It sets the two 19F shift entries to `0.0000` and comments that the actual value is `-184.1865`, but is zeroed because that value does not matter here and the calculation is faster. Scalar couplings are entered explicitly in the example; their units are not labelled.

The basis is `sphten-liouv` with `IK-0` and `inter_level=1`. Three manual fragment memberships are specified, with `S3` symmetry for spin groups `[14 15 16]` and `[21 22 23]`; the basis also sets `longitudinal={{'19F'}}` and `projections={1}`. The code leaves ZTE off by default. A GPU enable statement is present only as a comment and is not active. No relaxation model is configured in this example.

The acquisition selects `{'1H'}`, sets the initial state to `state(spin_system,'L+','1H')`, and leaves the decoupling list empty. Offset is `1400`, sweep `2500`, acquired points `4096`, zero fill `16536`, and the axis is in ppm with `invert_axis=1`. The code does not label units for the field, offset, or sweep literals. The receiver uses the same operator description with `coil_state` instead.

## Propagation and processing

The pulse-acquire sequence is propagated with `liquid(...,@acquire,...,'nmr')`. The FID is apodised with `exp` and parameter `5`, Fourier-transformed, and plotted as the real spectrum with the frequency axis inverted. The source supplies no numerical peak positions or intensities.

## Source

[MATLAB example](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_liquids/pa_difluoroheptane_anti.m)
