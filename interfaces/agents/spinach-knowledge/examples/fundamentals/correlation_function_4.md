# examples/fundamentals/correlation_function_4.m

- MATLAB implementation: [examples/fundamentals/correlation_function_4.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/fundamentals/correlation_function_4.m)

- Signature: `correlation_function_4()`

## What is compared

This high-rank isotropic case compares a Monte Carlo estimate with Spinach's analytical expansion for `G(L,k,m,p,q)=<D{L}(k,m)*D{L}(p,q)'>`. It sets `sigma_iso=0.2`, `L=4`, and `k=-1, m=2, p=-1, q=2`; each Wigner index is shifted by `L+1` for MATLAB indexing.

## Numerical construction

The source generates `1e6` direction-cosine-matrix updates from three independent standard-normal increments, with all three displayed skew-symmetric generators scaled by `sigma_iso`. It stores the trajectory, converts each matrix through `dcm2euler` to `wigner(L,...)` in a `parfor` loop, and computes normalised `xcorr` for the selected elements. The result is rescaled by `1/(2*L+1)` (here `1/9`) and shifted together with its lag vector by `ifftshift`. With `nlags=100`, the plotted portion corresponds to lag points 0 through 99.

## Analytical construction and observable

A one-spin dummy system (zero magnet, isotope `G`) is configured with Redfield relaxation, lab-frame retention, zero equilibrium, spherical-tensor Liouville formalism, and no basis approximation. It uses `tau_c={1/(3*sigma_iso^2)}`. The output of `corrfun(spin_system,L,k,m,p,q)` is summed as `sum_j weights{1}(j)*exp(rates{1}(j)*(0:(nlags-1)))`. The function plots real Monte Carlo points against the analytical curve; it does not set a tolerance or assert pass/fail. The source estimates minutes of calculation time.
