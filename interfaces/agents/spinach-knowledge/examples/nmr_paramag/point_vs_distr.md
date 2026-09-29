# examples/nmr_paramag/point_vs_distr.m

- Signature: `point_vs_distr()`
- Source: [examples/nmr_paramag/point_vs_distr.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_paramag/point_vs_distr.m)
- Model reference cited by the source: [DOI 10.1039/c6cp05437d](https://doi.org/10.1039/c6cp05437d)

## Synthetic density and PCS data

The example constructs four equal-weight Gaussian components with randomly selected centroids inside a 3-Angstrom cube and `sigma=0.5 Angstrom`. It places 100 nuclei randomly in a spherical layer with inner radius 5 Angstrom and thickness 10 Angstrom; the grid has 64 points per axis over the corresponding -15 to 15 Angstrom extent. The susceptibility tensor is rotated by Euler angles (pi/3, pi/4, pi/5) from the diagonal values formed with `ax=-0.45` and `rh=-0.05`; the source labels these tensor values in cubic Angstroms.

After padding the density with two original grid lengths on each side, `kpcs(...,'fft')` generates the PCS values at those nuclei. `ippcs` fits the point model about the origin; `ilpcs` fits multipole ranks 0, 1, and 2 about the same origin. The source prints the true tensor, fitted parameters and their reported standard deviations, plots residuals against distance, and compares true and fitted multipole moments.

## Scope and omissions

This is randomised synthetic data, with no fixed seed, experimental spectrum, field, or temperature. It is not a carbonic-anhydrase case and has no named metal site or protein residue. The source reports outputs at run time rather than fixed fit values; its header estimates a runtime of minutes.
