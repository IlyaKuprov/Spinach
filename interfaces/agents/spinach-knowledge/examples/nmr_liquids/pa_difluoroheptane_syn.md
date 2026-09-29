# examples/nmr_liquids/pa_difluoroheptane_syn.m

- Signature: `pa_difluoroheptane_syn()`
- Source: [MATLAB example](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_liquids/pa_difluoroheptane_syn.m)

## Purpose

Simulates a liquid-state, pulse-acquire 1H NMR FID for syn-3,5-difluoroheptane. The example is not an INADEQUATE, inversion-recovery, NOE, or NOESY experiment: it prepares and detects 1H single-quantum magnetisation with `liquid(...,@acquire,...,'nmr')`. The molecule's 19F spins are included in the model, but the acquisition channel, initial state, and receiver are all 1H.

## System and basis

The source defines 23 spins: seven 12C, fourteen 1H, and two 19F, at `sys.magnet=11.7464` T. It specifies scalar shifts and scalar couplings in the source; examples of the coded shift values are 1.0189, 4.6138, and 0.0000. The basis is manually partitioned into three fragment subspaces (`bas.manual`, `inter_level=1`) in the `sphten-liouv` formalism with `IK-0` approximation. It applies S3 symmetry to the two three-proton groups, includes 19F longitudinal terms, and sets projection 1. The source disables ZTE; its GPU enable line is commented out.

## Acquisition and processing

The source sets `parameters.spins={'1H'}`, with 1H `L+` initial state and receiver, and no decoupling. Acquisition settings are `offset=1400`, `sweep=2500`, `npoints=4096`, and `zerofill=16536`; the plotted axis is ppm and is inverted. No units for offset or sweep are stated in this source. It applies exponential apodisation with parameter 5, computes the zero-filled shifted Fourier transform, takes its real part, and plots it.

## Source limits

The source cites [the reported synthesis and NMR study](https://doi.org/10.1021/acs.joc.4c00670) and estimates minutes of calculation, faster with a GPU. It provides a spin model and processing recipe, not a measured spectrum, reported rate, or peak list. The source code does not specify an explicit relaxation model.
