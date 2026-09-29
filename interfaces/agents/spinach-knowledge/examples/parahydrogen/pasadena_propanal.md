# examples/parahydrogen/pasadena_propanal.m

- Signature: `pasadena_propanal()`
- Source: [examples/parahydrogen/pasadena_propanal.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/parahydrogen/pasadena_propanal.m)

## Purpose and PASADENA context

The source describes a PASADENA simulation of parahydrogenation of acrolein into propanal; it gives a calculation time of seconds and credits Ronghui Zhou and Ilya Kuprov. Unlike the ethylbenzene source, this file gives no DOI. In the physical experiment, parahydrogen's proton singlet order becomes observable after chemical addition makes the product proton environments inequivalent. This script does not simulate that reaction or the singlet-to-product transfer: its initial state is directly set to longitudinal two-spin order on product spins 1 and 4. It is not an ALTADENA low-field-to-high-field transfer or SABRE catalyst-exchange simulation.

## Product spin system

The six-spin model contains only `1H` nuclei at `7.05` T. Under Spinach's nuclear Zeeman and scalar-coupling unit conventions, the shifts in index order are `{1.11, 1.11, 1.11, 2.46, 2.46, 9.79}` ppm. The assigned scalar couplings are `J(1-4)=J(2-4)=J(3-4)=J(1-5)=J(2-5)=J(3-5)=7.3` Hz and `J(4-6)=J(5-6)=1.4` Hz; other matrix entries are not explicitly assigned. The basis is spherical-tensor Liouville space with no basis approximation and symmetry groups `S3` on spins 1–3 and `S2` on spins 4–5. The initial state is `state(spin_system,{'Lz','Lz'},{1,4})`.

## Acquisition and observable

The proton acquisition uses a `pi/4` y pulse, 500 Hz transmitter offset, 1000 Hz sweep, and 1024 points. The FID is zero-filled to 8192 points, Gaussian-apodised with parameter 10, Fourier transformed, and plotted as a real spectrum with the ppm axis inverted. No catalyst, chemical-exchange kinetics, explicit parahydrogen singlet, or relaxation superoperator appears in the script. This is a calculated product spectrum, not a measured trace or a prediction of conversion or polarisation yield.
