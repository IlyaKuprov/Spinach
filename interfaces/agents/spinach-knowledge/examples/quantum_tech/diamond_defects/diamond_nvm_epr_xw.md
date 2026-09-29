# examples/quantum_tech/diamond_defects/diamond_nvm_epr_xw.m

[Source MATLAB example](https://github.com/IlyaKuprov/Spinach/blob/main/examples/quantum_tech/diamond_defects/diamond_nvm_epr_xw.m)

## Physical model

This is a conventional powder EPR field sweep for the negatively charged nitrogen-vacancy centre in diamond, using its ground-state model. The example chooses 14N and orientation 111. The model builder defines an electron spin triplet (Spinach electron label `E3`, S = 1) coupled to a 14N nucleus (I = 1); its spin Hamiltonian contains electron zero-field splitting, nitrogen hyperfine coupling, and nitrogen quadrupolar coupling. `parameters.spins={'E3'}` selects the electron-spin EPR transitions, with nuclear-spin structure affecting their positions and intensities.

No cavity mode, Jaynes–Cummings/Tavis–Cummings interaction, vacuum-Rabi coupling, detuning, driven cavity dynamics, or dissipation appears in this example. The microwave frequency is used for the resonance condition in a magnetic-field sweep; the script specifies no drive amplitude or time-dependent observable. This is not a cavity-device simulation and does not imply measured performance or device fidelity.

## Calculation and plotted quantity

The full Zeeman Hilbert basis is used (`zeeman-hilb`, `approximation='none'`), with spherical grid `rep_2ang_100pts_sph`. The source sets `sys.magnet` to 1 T, and `fieldsweep` runs at 9.5 GHz over 0.1–0.5 T and at 94 GHz over 3.2–3.5 T. Both sweeps use 0.001 T FWHM, integration tolerance 1e-4, transition-moment tolerance 0.01, RSPT order `Inf`, and 512 field points. The plotted arrays are simulated EPR intensity in arbitrary units against returned magnetic-field axes, not experimental spectra or time traces. This example plots the `fieldsweep` output rather than declaring a separate observable operator.

## Parameter provenance

The ground-state NV magnetic parameters are attributed by `diamond_nvm_gs` to Felton et al., *Physical Review B* 79, 075203 (2009), [doi:10.1103/PhysRevB.79.075203](https://doi.org/10.1103/PhysRevB.79.075203).