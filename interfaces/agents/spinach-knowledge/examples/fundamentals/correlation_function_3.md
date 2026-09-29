# examples/fundamentals/correlation_function_3.m

- MATLAB implementation: [examples/fundamentals/correlation_function_3.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/fundamentals/correlation_function_3.m)

- Signature: `correlation_function_3()`

## What is compared

This example compares a finite-trajectory estimate with Spinach's analytical correlation-function expansion for `G(L,k,m,p,q)=<D{L}(k,m)*D{L}(p,q)'>`. It tests a rhombic rotational-diffusion model at `L=2`, with `k=-1, m=2, p=-1, q=2` and `sigma_x=0.1`, `sigma_y=0.2`, `sigma_z=0.3`. As in the source, each index is shifted by `L+1` before it indexes a Wigner matrix.

## Numerical construction

The trajectory contains `1e6` updates and the correlation uses `nlags=300`. Three standard-normal increments weight the three explicit skew-symmetric generators; the source orders those terms as `sigma_z` times the first matrix, `sigma_y` times the second, and `sigma_x` times the third. The combined increment is exponentiated and right-multiplied into a direction-cosine matrix initialised to the identity. A `parfor` loop converts each stored matrix by `dcm2euler` and `wigner(L,...)`; normalised `xcorr` of the selected Wigner elements is multiplied by `1/(2*L+1)`. `ifftshift` places zero lag first, and the plotted segment uses lag points 0 through 299.

## Analytical construction and observable

The dummy spin system has zero magnet, isotope `G`, Redfield relaxation, lab-frame retention, zero equilibrium, spherical-tensor Liouville formalism, and no basis approximation. Its correlation-time array is set to `1./(3*[sigma_x sigma_y sigma_z].^2)` (x/y/z order), whereas the Monte Carlo generator terms above are written z/y/x; keep these source orderings distinct. The weights and rates from `corrfun(spin_system,L,k,m,p,q)` form `sum_j weights{1}(j)*exp(rates{1}(j)*(0:(nlags-1)))`. The function plots the real Monte Carlo estimate against the analytical curve and contains no tolerance or pass/fail assertion. The source estimates minutes of calculation time.
