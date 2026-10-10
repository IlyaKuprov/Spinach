# examples/nmr_solids/case_studies/mathies_carbonate/sle_nmr_dd_csa_mhc.m

Source: [examples/nmr_solids/case_studies/mathies_carbonate/sle_nmr_dd_csa_mhc.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_solids/case_studies/mathies_carbonate/sle_nmr_dd_csa_mhc.m)

## Purpose and spin model

This example calculates a proton spectrum for two water protons in a monohydrocalcite unit cell, including isotropic rotational diffusion and magic-angle spinning (MAS). It reads CASTEP data from
yt `mhc.magres`, removes C, O, and Ca, and keeps the proton sites at entries 1 and 4. Their shielding tensors are converted to shift tensors as `29.25*eye(3)-props.cst{1}` and `29.25*eye(3)-props.cst{4}`; the corresponding coordinates are taken from those same entries. The source attributes this shift parametrisation to Huang et al., ACIE 2021, and cites [10.1038/s41467-023-44381-x](https://doi.org/10.1038/s41467-023-44381-x).

The source labels the experiment 400 MHz NMR and sets `sys.magnet=9.4`. It uses the full `sphten-liouv` basis with no approximation. The MAS parameters are `rate=10000` and `axis=[1 1 1]`; the source does not attach a unit to the rate value. The initial state and receive coil are both proton `L+`, and acquisition is configured for `1H` with no decoupled spins.

## Acquisition and displayed result

The simulation sweeps five rotational-diffusion correlation times, `1e-7` to `1e-3` s, with Wigner ranks `[2 3 5 7 13]`, respectively. For each value it calls `gridfree` with `@acquire` and the NMR mode, using 1024 points, a sweep value of 120000, zero fill to 4096, and offset 0. The source does not state units for those three acquisition values. It applies exponential apodisation with parameter 6, Fourier transforms the FID, and plots the real spectrum with amplitude in arbitrary units and a displayed vertical range from -1 to 10. The legend identifies the five correlation times.

This is a simulation, not a spectrum read from an experimental acquisition: the input is the CASTEP file and the programmed spin/acquisition model. The output is the calculated spectrum series. The source estimates minutes of calculation, or seconds with a GPU; GPU enablement is commented out in this file. It does not provide numerical peak positions or a measured-versus-simulated comparison.

Zero track elimination is explicitly enabled with `zte` in `sys.enable`.
