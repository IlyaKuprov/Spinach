# examples/nmr_solids/mas_powder_gly_floquet.m

Source: [examples/nmr_solids/mas_powder_gly_floquet.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_solids/mas_powder_gly_floquet.m)

## Model

This is a simulated 13C MAS powder example for glycine. It reads `../standard_systems/glycine.log` with `gparse` and passes the parsed data to `g2spinach`, selecting 13C and 15N spins with reference values [182.1, 264.5]. The code labels the input as PCM DFT-derived spin-system properties and sets the field to 14.1 T. The source does not state units for the two reference values or for the interaction and proximity cutoffs (5.0 and 4.0). These are model-building inputs, not measured spectrum values.

The basis is spherical-tensor Liouville space with no approximation, projection +1, and a longitudinal 15N specification. The experiment settings are a [1, 1, 1] rotor axis, rate 2000 Hz, maximum rank 23, and grid `leb_2ang_rank_23`. The header says the spectrum assumes proton decoupling, but the selected spin list is 13C and 15N and `parameters.decouple` is empty; the source does not configure an explicit proton spin or RF decoupling field. The simulation call is `floquet(spin_system,@acquire,parameters,'nmr')`.

## Signal and display

Both the initial state and receiver coil are the 13C `L+` state. The code requests 256 points, a sweep setting of 5e4, zero filling to 1024, offset 17000, and a ppm axis with inversion enabled. It applies exponential apodisation parameter 6 to the FID, Fourier transforms with `fftshift`, then plots the real spectrum using `plot_1d`. The result is a computed spectrum from the imported model and acquisition settings, not an experimentally measured spectrum.
