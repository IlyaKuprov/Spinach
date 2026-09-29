# examples/quantum_tech/diamond_defects/diamond_siv0_epr_xw.m

[Example source](https://github.com/IlyaKuprov/Spinach/blob/main/examples/quantum_tech/diamond_defects/diamond_siv0_epr_xw.m) · [SiV0 spin-system builder](https://github.com/IlyaKuprov/Spinach/blob/main/etc/diamond_defects/diamond_siv0.m)

- Signature: `diamond_siv0_epr_xw()`

## Model

This example calculates powder EPR field sweeps for the neutral silicon-vacancy centre in diamond. The call to `diamond_siv0` selects `silicon='29Si'`, [111] orientation, and `n_13c=0`: its model therefore includes the SiV0 spin-1 electronic state and one ²⁹Si nuclear spin, but no added ¹³C nuclei. The helper's electron g principal values are 2.0035, 2.0035 and 2.0042; it sets a zero-field-splitting parameter of 1.000 GHz and a ²⁹Si hyperfine tensor with principal components 78.9, 78.9 and 76.3 MHz. The EPR routine is called with `parameters.spins={'E3'}`, selecting the spin-1 electron transition channel while retaining the silicon hyperfine interaction in the Hamiltonian.

The magnetic parameters in the helper are attributed to Edmonds et al., *Phys. Rev. B* 77, 245205 (2008), [doi:10.1103/PhysRevB.77.245205](https://doi.org/10.1103/PhysRevB.77.245205).

## Field-swept spectra

The full Zeeman-Hilbert basis is used without an approximation. The powder grid is `rep_2ang_100pts_sph`. The X-band scan uses 9.5 GHz and 0.1–0.5 T; the W-band scan uses 94 GHz and 3.2–3.5 T. Each spectrum is calculated at 2,048 field points and plotted as simulated intensity in arbitrary units versus magnetic field in tesla. These are model spectra, not measured defect data.

The source sets `fwhm=0.001` but does not annotate the unit, so no unit conversion is asserted. The page does not claim a measured line shape, runtime benchmark, or convergence result.
