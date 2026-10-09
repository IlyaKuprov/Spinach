# examples/parahydrogen/pasadena_ethylbenzene.m

- Signature: `pasadena_ethylbenzene()`
- Source: [examples/parahydrogen/pasadena_ethylbenzene.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/parahydrogen/pasadena_ethylbenzene.m)

## Purpose and PASADENA context

The source describes a PASADENA simulation for parahydrogenation of styrene to ethylbenzene, intended to reproduce the top trace of Figure 5 in [doi:10.1039/b914188j](https://doi.org/10.1039/b914188j). It gives a calculation time of seconds. In the physical experiment, the two protons from parahydrogen begin in a nuclear-spin singlet; chemical addition to an unsaturated precursor makes their product environments inequivalent, allowing the non-equilibrium spin order to appear in the NMR signal. Here the reaction and singlet-to-product transfer are not propagated: the code begins directly with a chosen product-spin density operator, `Lz` order on spins 1 and 4. It is not an ALTADENA low-field-to-high-field transfer or SABRE catalyst-exchange simulation.

## Product spin system

The model has ten `1H` spins at `7.05` T. Under Spinach's nuclear Zeeman and scalar-coupling unit conventions, the shifts in spin-index order are `{1.201, 1.201, 1.201, 2.625, 2.625, 7.207, 7.207, 7.265, 7.265, 7.155}` ppm. The scalar-coupling entries assigned in the source are `J(1-4)=J(2-4)=J(3-4)=J(1-5)=J(2-5)=J(3-5)=7.63` Hz; `J(6-8)=7.63` Hz; `J(6-10)=J(7-10)=1.26` Hz; `J(7-9)=7.63` Hz; and `J(8-10)=J(9-10)=7.44` Hz. These are the nonzero assignments written in the source; it does not provide an atom-numbering diagram here.

The basis uses spherical-tensor Liouville space, the `IK-2` approximation, scalar-coupling connectivity with proximity level 1, and permutation groups `S3` on spins 1–3 and `S2` on spins 4–5. The initial state is `state(spin_system,{'Lz','Lz'},{1,4})`.

## Acquisition and observable

The code simulates a proton FID with a `-pi/4` y pulse, 500 Hz transmitter offset, 1000 Hz sweep, and 1024 acquired points. It zero-fills to 8192 points, applies Gaussian apodisation with parameter 10, Fourier transforms, and plots the real spectrum with the ppm axis inverted. This is a computed PASADENA-style product spectrum; no catalyst, hydrogenation time course, exchange rate, or relaxation superoperator is specified in the script, and the plotted trace is not itself an experimental measurement.
