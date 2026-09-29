# examples/fitting/fluoroalkanes/anti_difluoroheptane.m

- MATLAB implementation: [examples/fitting/fluoroalkanes/anti_difluoroheptane.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/fitting/fluoroalkanes/anti_difluoroheptane.m)

- Signature: `anti_difluoroheptane()`

## Purpose

Fit one 19F and two 1H spectral regions for anti-3,5-difluoroheptane by varying scalar couplings and spectrum scales. This is broader than the source header's 1H-only description: the implementation explicitly loads and fits the fluorine spectrum as well. The paper cited by the source is [here](https://doi.org/10.1021/acs.joc.4c00670); the source estimates hours of calculation time.

## Data and spin model

The function loads `anti_dfh_fluorine.mat`, `anti_dfh_proton_a.mat`, and `anti_dfh_proton_b.mat`; each supplies `axis_ppm` and `spec`. Each experimental spectrum is divided by its own maximum. The 23-site model at field `11.7464` includes 1H, 19F, 12C, and `G` ghost sites. The source comment says the methyl groups are ghosted because they do not influence the selected signals; the code assigns fixed couplings of `7.45` to the ghosted methyl sites.

## Parameterisation and objective

The starting vector is `[-14.4404 2.1802 10.0564 14.0478 38.0147 27.7945 18.4317 1.2295 3.1159 -15.1172 11.2721 49.6307 18.9158 4.5308 7.5739]`. Twelve entries set selected scalar couplings and three scale the 19F and two 1H signals. The source uses `fminsearch` with `MaxIter=5000` and `MaxFunEvals=Inf`. Its objective is the squared 19F residual norm plus twice the squared residual norm for each 1H spectrum.

## Acquisition and display

The function makes separate, non-decoupled liquid acquisitions for 19F and the two 1H signal selections. All three declare ppm axes with reversed direction. The 19F sequence sets offset `-86700`, sweep `300`, 512 acquired points, and zero filling to 2048; the 1H sequences set offsets `2335` and `890`, sweep `128` each, 256 acquired points, and zero filling to 1024. The source does not label the numeric offset units. The 19F FID receives exponential apodisation, while the 1H FIDs receive Gaussian apodisation. After Fourier transformation, the simulated axes are built from the sequence settings and the spectra are interpolated to the experimental axes with `pchip`. Each objective evaluation plots experiment and simulation in three stacked panels; the function also displays the parameter vector and final answer.

This is a source-described local search, not a report of a measured optimum: the entry point has no output argument and the source supplies no uncertainty estimate or fit-quality result. The 1H residuals receive twice the weight of the 19F residual in the objective, so the three plotted panels do not contribute equally to the stated error functional.
