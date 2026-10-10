# examples/esr_sol_swept/temperature_gd_dota.m

- Signature: `temperature_gd_dota()`
- Source: [`examples/esr_sol_swept/temperature_gd_dota.m`](https://github.com/IlyaKuprov/Spinach/blob/main/examples/esr_sol_swept/temperature_gd_dota.m)

## Purpose

A compact temperature series of powder-averaged, W-band, field-swept ESR spectra for the Gd(III)-DOTA model. It is useful as a worked example of setting temperature while calculating powder ESR spectra for a one-electron-spin model.

## Spin system and model

The source declares one `E8` spin, sets the isotropic Zeeman value to 1.9918, and supplies a zero-field-splitting tensor through `inter.coupling.eigs{1,1}=[0.57e9 0.57e9 -2*0.57e9]/3` with Euler angles `[0 0 0]`. The example identifies this as Gd(III) DOTA. It uses the `zeeman-hilb` formalism with `bas.approximation='none'`, i.e. the script's exact, unapproximated Hilbert-space basis choice.

## Sweep and temperature settings

The four temperatures are `[100 10 1 0.1]` K (the source labels these values Kelvin). At each point the script assigns `inter.temperature`, rebuilds the spin system and basis, then calls `fieldsweep` with spin `E8`, grid `rep_2ang_100pts_sph`, `mw_freq=90e9`, `fwhm=2e-4`, `int_tol=1.0`, `tm_tol=0.1`, field window `[3.05 3.4]`, `npoints=4096`, and `rspt_order=Inf`. The script does not annotate units for the numeric microwave-frequency, linewidth, tolerance, or resonance-order fields; it does explicitly label the plotted field axis in tesla.

There is no pulse-program sequence in this example: the observable is calculated by the field-sweep ESR engine, not by time-domain pulse evolution. A 2-by-2 figure places one spectrum at each temperature, with magnetic field in tesla and intensity in arbitrary units.

## Output and scope

The result is a comparison of four simulated powder spectra, not a reported experimental fit or a tabulated peak assignment. The source comment estimates calculation time as seconds; that is the example's note, not a benchmark guarantee. The code contains no hyperfine couplings and no DOI or literature citation.
