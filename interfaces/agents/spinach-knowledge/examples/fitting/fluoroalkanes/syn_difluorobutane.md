# examples/fitting/fluoroalkanes/syn_difluorobutane.m

- Signature: `syn_difluorobutane()`

## Purpose

Fit the 1H NMR spectrum of syn-2,3-difluorobutane by varying J-couplings, Gaussian linewidth, and spectral amplitude. For further details, see https://doi.org/doi/10.1021/acs.joc.4c00670. Calculation time: hours.

## Workflow

- Load the CH and ME experimental frequency axes and spectra from `syn_dfb_proton.mat`. Normalize the two spectra by their respective integrals to −2 and −6, then concatenate the intervals for fitting.
- Start from the parameter vector `[23.95 6.47 0.90 4.36 18.15 47.88 -11.61 13.63 1.7]`, ordered as near F–CH3, CH3–H, far F–CH3, three-bond H–H, three-bond F–H, two-bond F–H, three-bond F–F couplings, Gaussian linewidth, and amplitude.
- Use `fminsearch` to minimize the squared norm of the difference between the real experimental and simulated spectra. Optimization uses central finite differences, `DiffMinChange=1e-3`, `MaxIter=5000`, and `MaxFunEvals=Inf`; the fitted parameters are displayed.

## Spectrum simulation and comparison

- Construct an eight-1H/two-19F spin system at 11.7464 T. Assign chemical shifts of 1.34375 to the six methyl protons, 4.57625 to the two methine protons, and 0.00 to both fluorines. The seven coupling parameters specify the symmetry-related scalar couplings between methyl protons, methine protons, and fluorines.
- Use an untruncated Zeeman–Hilbert basis with separate S3 symmetry groups for the two methyl groups. Acquire a 1H spectrum with `L+` initial and detection states, no decoupling, an offset of 1500 Hz, an 1800 Hz sweep, and 4096 points.
- Apply Gaussian apodisation using the fitted linewidth, Fourier-transform with zero filling to 32768 points, scale by the fitted amplitude, reverse the spectrum, and interpolate it onto the concatenated experimental frequency axes using `pchip`.
- Plot experimental and simulated real spectra in the 2230–2350 Hz and 645–700 Hz windows during optimization.