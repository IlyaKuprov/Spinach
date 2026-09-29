# examples/quantum_tech/diamond_defects/diamond_p1_epr_xw.m

[Example source](https://github.com/IlyaKuprov/Spinach/blob/main/examples/quantum_tech/diamond_defects/diamond_p1_epr_xw.m) · [P1 spin-system builder](https://github.com/IlyaKuprov/Spinach/blob/main/etc/diamond_defects/diamond_p1.m)

- Signature: `diamond_p1_epr_xw()`

## Model

This example calculates field-swept powder EPR spectra for the diamond P1 centre, a substitutional nitrogen electron-spin defect. It calls `diamond_p1` with `nitrogen='14N'` and the [111] orientation. The model contains an electron spin S = 1/2 and the selected ¹⁴N nucleus (I = 1); the builder supplies the anisotropic electron g tensor, electron-nuclear hyperfine tensor and the ¹⁴N quadrupolar interaction. For this isotope, the helper gives principal g values 2.00220, 2.00220 and 2.00218, hyperfine components 81.3, 81.3 and 114.0 MHz, and a quadrupolar parameter of −3.97 MHz. The spin-system parameters cite Nir-Arad et al., *Phys. Chem. Chem. Phys.* 26, 27633 (2024), [doi:10.1039/d4cp03055a](https://doi.org/10.1039/d4cp03055a), and Smith et al., *Phys. Rev.* 115, 1546 (1959), [doi:10.1103/PhysRev.115.1546](https://doi.org/10.1103/PhysRev.115.1546).

The example uses the exact Zeeman-Hilbert-space basis (`zeeman-hilb`, no approximation) and asks the EPR field-sweep routine for electron transitions with `parameters.spins={'E'}`. That selects the electron-spin transitions; the nitrogen coupling shapes their hyperfine structure rather than being selected as an independent EPR-active spin.

## Field-swept spectra

The powder average uses `rep_2ang_100pts_sph`. The X-band calculation uses 9.5 GHz microwave frequency and scans 0.330–0.350 T with 1,024 field points. The W-band calculation uses 94 GHz and scans 3.348–3.360 T, reusing the same point count. The two plotted traces are simulated intensity in arbitrary units versus magnetic field in tesla. They are model spectra, not measured defect data.

The example sets `fwhm=1e-4`, but does not annotate that input's unit; it is therefore not converted here. No measured spectrum, runtime benchmark, or convergence study is reported.
