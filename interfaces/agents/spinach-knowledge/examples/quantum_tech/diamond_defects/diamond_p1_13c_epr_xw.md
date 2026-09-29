# examples/quantum_tech/diamond_defects/diamond_p1_13c_epr_xw.m

[Source MATLAB example](https://github.com/IlyaKuprov/Spinach/blob/main/examples/quantum_tech/diamond_defects/diamond_p1_13c_epr_xw.m)

## Physical model

This is a conventional powder EPR field sweep for the P1 substitutional-nitrogen centre in 13C-enriched diamond. The example selects 14N and orientation 111. Its spin system contains an electron doublet (electron label `E`, S = 1/2), the 14N nucleus (I = 1), and neighbouring 13C nuclei (I = 1/2); the 14N has hyperfine and quadrupolar interactions. The builder starts with 18 13C spins, and the example retains the nitrogen plus only 13C spins whose isotropic hyperfine couplings exceed 8 MHz. The EPR transition family is selected with `parameters.spins={'E'}`; nitrogen and retained carbon hyperfine couplings shape the spectrum.

There is no cavity mode, Jaynes–Cummings/Tavis–Cummings or vacuum-Rabi model, cavity coupling, detuning, driven coherent evolution, or dissipation in this source. The microwave frequency sets the EPR resonance condition for a field sweep; no drive amplitude is specified. The calculation is not a device-fidelity estimate or a claim about measured defect spectra.

## Calculation and plotted quantity

The full Zeeman Hilbert basis is used (`zeeman-hilb`, `approximation='none'`), with spherical grid `rep_2ang_100pts_sph`. The source sets `sys.magnet` to 1 T and uses 9.5 GHz over 0.31–0.36 T and 94 GHz over 3.33–3.38 T. Both sweeps use 5e-4 T FWHM, integration tolerance 0.1, transition-moment tolerance 0.1, RSPT order `Inf`, and 1024 field points. The plotted output is simulated EPR intensity in arbitrary units versus the returned magnetic-field axes, rather than a time-domain signal or measurement. The script does not declare an explicit observable matrix; it plots the spectra returned by `fieldsweep`.

## Parameter provenance

The electron and nitrogen parameters are inherited from `diamond_p1`, which cites Nir-Arad et al. (2024), [doi:10.1039/d4cp03055a](https://doi.org/10.1039/d4cp03055a), and Smith et al. (1959), [doi:10.1103/PhysRev.115.1546](https://doi.org/10.1103/PhysRev.115.1546). The 13C hyperfine tensors and site assignments in `diamond_p1_13c` are attributed to Barklie and Guven (1981), [doi:10.1088/0022-3719/14/25/009](https://doi.org/10.1088/0022-3719/14/25/009); Cox, Newton, and Baker (1994), [doi:10.1088/0953-8984/6/2/012](https://doi.org/10.1088/0953-8984/6/2/012); and Peaker et al. (2016), [doi:10.1016/j.diamond.2016.10.013](https://doi.org/10.1016/j.diamond.2016.10.013).