# examples/quantum_tech/diamond_defects/diamond_ninter_epr_xw.m

[Source MATLAB example](https://github.com/IlyaKuprov/Spinach/blob/main/examples/quantum_tech/diamond_defects/diamond_ninter_epr_xw.m)

## Physical model

This is a conventional spin-Hamiltonian powder EPR simulation of the WAR9 nitrogen-interstitial centre in diamond. The example selects 15N, orientation 111, and the builder specifies one electron doublet (electron label `E`, S = 1/2) hyperfine-coupled to that nitrogen nucleus (I = 1/2). The WAR9 15N hyperfine principal values in `diamond_n_inter` are 8.30, 7.85, and 8.17 MHz, and its electron g tensor is anisotropic. The requested transition family is electron-spin EPR (`parameters.spins={'E'}`); the nuclear coupling contributes hyperfine structure.

This source has no cavity mode or Jaynes–Cummings/Tavis–Cummings or vacuum-Rabi term, cavity coupling, detuning, drive-amplitude dynamics, or dissipation. Its microwave frequency specifies the EPR resonance condition in a field sweep, not a coherent cavity drive. It therefore does not simulate a cavity device or establish measured spectra or device fidelity.

## Calculation and plotted quantity

The calculation uses the full Zeeman Hilbert basis (`zeeman-hilb`, `approximation='none'`) and spherical grid `rep_2ang_100pts_sph`. The source sets `sys.magnet` to 1 T and calls `fieldsweep` at 9.755 GHz over 0.347–0.349 T and at 94 GHz over 3.351–3.355 T. Both sweeps use 1e-5 T FWHM, integration tolerance 0.1, transition-moment tolerance 0.1, RSPT order `Inf`, and 1024 field points. The plotted observable is the simulated EPR intensity (arbitrary units) versus each returned field axis, not a measured trace or a time-domain response. No explicit observable matrix is constructed in this example; the plotted arrays come from `fieldsweep`.

## Parameter provenance

The WAR9 magnetic parameters are attributed by `diamond_n_inter` to Felton et al., *J. Phys.: Condens. Matter* 21, 364212 (2009), [doi:10.1088/0953-8984/21/36/364212](https://doi.org/10.1088/0953-8984/21/36/364212).