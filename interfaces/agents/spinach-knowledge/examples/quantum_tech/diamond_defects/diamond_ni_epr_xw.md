# examples/quantum_tech/diamond_defects/diamond_ni_epr_xw.m

[Source MATLAB example](https://github.com/IlyaKuprov/Spinach/blob/main/examples/quantum_tech/diamond_defects/diamond_ni_epr_xw.m)

## Physical model

This is a conventional spin-Hamiltonian, field-swept powder EPR calculation for the diamond Ni NE1 centre, not a cavity-QED simulation. The model builder identifies an electron doublet (Spinach electron label `E`, S = 1/2) coupled to two 14N nuclei; its NE1 parameter set has an anisotropic electron g tensor and anisotropic nitrogen hyperfine tensors. The example fixes the centre orientation as 111 before averaging over the spherical powder grid. EPR selection is on electron-spin transitions, requested with `parameters.spins={'E'}`; nuclear hyperfine structure shifts and splits those transitions.

There is no cavity mode, Jaynes–Cummings or Tavis–Cummings interaction, vacuum-Rabi splitting, cavity coupling, detuning parameter, driven time evolution, or dissipation model in this source. The microwave frequency is used to define the EPR resonance condition, not as a specified drive amplitude. Consequently these calculations say nothing about cavity coupling or device fidelity.

## Calculation and plotted quantity

The full Zeeman Hilbert basis is used (`zeeman-hilb`, `approximation='none'`). The source sets `sys.magnet` to 1 T, then `fieldsweep` calculates spectra over 0.30–0.36 T at 9.5 GHz and 3.1–3.4 T at 94 GHz. Both use `rep_2ang_100pts_sph`, a 1e-4 T FWHM, integration tolerance 1e-4, transition-moment tolerance 0.1, RSPT order `Inf`, and 2048 field points. The plotted arrays are simulated EPR intensity in arbitrary units against the returned magnetic-field axes; they are not experimental traces or time-domain signals. The example does not itself declare an operator matrix; it plots the spectrum returned by `fieldsweep`.

## Parameter provenance

The NE1 magnetic parameters are supplied by `diamond_ni`; its documentation cites Nadolinny et al., *Crystals* 7, 237 (2017), [doi:10.3390/cryst7080237](https://doi.org/10.3390/cryst7080237).