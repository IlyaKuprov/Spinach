# examples/parahydrogen/ortho_deuterium.m

- Signature: `ortho_deuterium()`

## Purpose

Simulates the ortho-deuteration spectrum of acrylonitrile shown in Figure 1 of the paper by Natterer, Greve, and Bargon ([doi:10.1016/S0009-2614(98)00784-2](https://doi.org/10.1016/S0009-2614(98)00784-2)). The source describes the calculation as taking seconds.

## Physical / mathematical content

The five-spin product contains three protons and two deuterons. The model uses the 4.697 T field labelled “Bargon's magnet” and scalar couplings specified for the deuterated product. In the Zeeman Hilbert-space basis, `deut_pair` constructs the deuteron singlet and five quintet states; their sum is used as the initial density operator for acquisition of the deuterium signal.

## Numerical / algorithmic content

A liquid-state acquisition applies a `pi/4` y pulse and simulates 1024 points over a 120 ppm sweep with an offset of 50 ppm. The FID is exponentially apodised (factor 6), zero-filled to 4096 points, Fourier transformed, and plotted.

## Implementation structure

The script defines the five isotopes and their scalar shifts and couplings, selects the Zeeman Hilbert formalism with no basis approximation, enables continuous deuteration, builds the singlet/quintet initial state, and sets the deuterium coil and pulse operators before acquisition and processing.
