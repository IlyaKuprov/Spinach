# examples/quantum_tech/diamond_defects/diamond_n2vm_epr_xw.m

- Source: [diamond_n2vm_epr_xw.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/quantum_tech/diamond_defects/diamond_n2vm_epr_xw.m)
- Signature: `diamond_n2vm_epr_xw()`

## Purpose and model

Calculate field-swept powder EPR spectra for the negatively charged N2V centre in diamond at X and W bands. The example configures `diamond_n2vm` with `orientation='111'`, `nitrogen='15N'`, and `include_13c=false`: both nitrogen nuclei use 15N and the optional 13C neighbours are omitted. It selects the electron EPR signal with `parameters.spins={'E'}`. The builder's magnetic parameters are attributed to Green et al., *Physical Review B* **92**, 165204 (2015), [doi:10.1103/PhysRevB.92.165204](https://doi.org/10.1103/PhysRevB.92.165204).

The calculation uses the `zeeman-hilb` basis without approximation and powder grid `rep_2ang_100pts_sph`. Its shared EPR inputs are `fwhm=0.00003`, `int_tol=0.1`, `tm_tol=0.1`, `npoints=1024`, and `rspt_order=Inf`.

## Sweep sequence and output

The source runs X band at 9.755 GHz (`mw_freq=9.755e9`) over 0.347–0.349 T, followed by W band at 94 GHz (`mw_freq=94e9`) over 3.351–3.355 T. Each returned `spec_x` or `spec_w` is plotted versus its field axis in tesla, with intensity in arbitrary units. These are simulated electron EPR spectra for the specified isotope model, not measurements of a diamond sample. The source comment estimates seconds of calculation time; this was not timed for this note.
