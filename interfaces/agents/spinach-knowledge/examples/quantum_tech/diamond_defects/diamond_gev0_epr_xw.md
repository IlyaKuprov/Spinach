# examples/quantum_tech/diamond_defects/diamond_gev0_epr_xw.m

- Source: [diamond_gev0_epr_xw.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/quantum_tech/diamond_defects/diamond_gev0_epr_xw.m)
- Signature: `diamond_gev0_epr_xw()`

## Purpose and model

Calculate field-swept powder EPR spectra for the neutral germanium-vacancy centre (GeV0) in diamond at X and W bands. The call to `diamond_gev0` sets `orientation='111'` and `germanium='none'`: the centre model retains its electron Zeeman and zero-field-splitting terms but does not add a germanium nuclear spin. The selected EPR signal is `parameters.spins={'E3'}`. Magnetic parameters for the builder are attributed to Nadolinny et al., *Phys. Status Solidi A* **213**, 2623 (2016), [doi:10.1002/pssa.201600211](https://doi.org/10.1002/pssa.201600211).

The simulation uses a `zeeman-hilb` basis without approximation and powder grid `rep_2ang_100pts_sph`. Shared settings are `fwhm=1e-3`, `int_tol=1e-3`, `tm_tol=0.01`, `npoints=2048`, and `rspt_order=Inf`; the source assigns these values and does not describe them as fitted experimental parameters.

## Sweep sequence and output

The model is swept at X band (9.5 GHz, `mw_freq=9.5e9`) over 0.05–0.45 T, then at W band (94 GHz, `mw_freq=94e9`) over 3.25–3.46 T. Each `fieldsweep` returns a spectrum and field axis; the plots use magnetic field in tesla and label intensity in arbitrary units. The outputs are calculated spectra, not recorded GeV0 measurements. The MATLAB source estimates seconds of calculation time; that estimate was not timed for this note.
