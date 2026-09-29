# examples/relaxation_theory/maz_noesy_2.m

[Source file](https://github.com/IlyaKuprov/Spinach/blob/main/examples/relaxation_theory/maz_noesy_2.m) · Signature: `maz_noesy_2()`

## Purpose and physical model

This example calculates a two-dimensional methylaziridine NOESY spectrum with two scalar-relaxation contributions: first-kind relaxation from nitrogen-centre inversion modulating scalar couplings, and second-kind relaxation associated with rapid quadrupolar relaxation of the `14N` nucleus. The source cites [the study describing the effect](https://doi.org/10.1002/ange.201410271) and gives an estimated calculation time of minutes.

The spin system is seven `1H` nuclei and one `14N` nucleus, with magnet setting 11.75. The source specifies vacuum-DFT shielding tensors, a vacuum-DFT quadrupole tensor for nitrogen, scalar couplings, and molecular coordinates; the coordinate comment explicitly gives Angstrom. The nitrogen quadrupole matrix is the source matrix [-1.2932, 0.6251, 1.8700; 0.6251, 1.7170, -2.3127; 1.8700, -2.3127, -0.4238] multiplied by 1e6, with no unit stated for its entries. Nitrogen-to-proton scalar-coupling entries include 4.4, 5.2, and 44.8, also without units stated. It assigns isotropic shifts from experiment using the values 1.3, 1.7, 1.9, 0.0, 0.1, 1.2, 1.2, and 1.2. This experimental assignment is an input to the calculation, not a measured spectrum or measured intensity supplied to the plotting step.

## Relaxation and sequence

The basis is `sphten-liouv` with IK-2 approximation, scalar-coupling connectivity, and proximity level 4. The source sets inter-spin and proximity cutoffs to 2.0 and 4.0 and disables Krylov propagation. It requests Redfield, SRFK, and SRSK relaxation, uses zero equilibrium, sets `srsk_sources` to 4, and keeps the relaxation terms in the kite. The Redfield correlation-time input is 25e-12; the SRFK inputs are `srfk_tau_c={[1.0 1e-3]}` and modulation depth 15.0 for entries (1,5), (2,5), and (3,5). The source does not label units for these relaxation parameters, so their supplied numeric values are retained without converting them to named units.

The source runs `liquid(spin_system,@noesy,parameters,'nmr')`. The `@noesy` sequence observes protons, starts from proton `Lz`, and uses a 2.0 s mixing time, 500 Hz offset, sweeps of 1400 in both dimensions, 256 points per dimension, and zero-filling to 1024 in each dimension. The cosine and sine acquisition components are cosine-apodised in both dimensions, Fourier transformed along F2, combined into a states signal, then Fourier transformed along F1. The plot is the negative real part of the calculated spectrum, shown with `plot_2d`. These lines and intensities are simulation output; the example does not load or plot a measured NOESY spectrum.
