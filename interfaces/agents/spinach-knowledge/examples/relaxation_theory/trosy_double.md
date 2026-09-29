# examples/relaxation_theory/trosy_double.m

- Signature: `trosy_double()`

## Purpose

A 13C-detected simulation of the source's Double TROSY example. The source estimates a calculation time of seconds. The four-spin fragment contains two `1H`, one `19F`, and one `13C`; its tensors, coordinates, and scalar-coupling submatrix are selected from a parsed DFT output. The DFT entries selected are `[10 20 19 8]` in that spin order. The source labels the input a 3-fluorotyrosine calculation, while the input filename is `4_fluoro_phe.out`; the structure's naming is therefore left as the source presents it rather than reconciled here.

## Relaxation and acquisition

The model uses secular Redfield relaxation, zero equilibrium, and `tau_c = 20e-9`. The field parameter is set to `14.1`; the source does not annotate its unit. The basis is the full `sphten-liouv` basis without approximation. The initial density operator and receiver are both the `13C` `L+` state, with no decoupled spins. The source sets offset `26800`, sweep `500`, 2048 points, zero filling to 16384 points, and a ppm axis with the displayed axis inverted.

The source acquires a liquid-state NMR signal, applies Gaussian apodisation with parameter `6`, Fourier transforms it, and plots the real spectrum. It then rebuilds the system after clearing both proton coordinates and repeats the acquisition; the source describes this comparison as removing proton dipolar relaxation. These are simulated spectra, not measurements.

Source: [examples/relaxation_theory/trosy_double.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/relaxation_theory/trosy_double.m).
