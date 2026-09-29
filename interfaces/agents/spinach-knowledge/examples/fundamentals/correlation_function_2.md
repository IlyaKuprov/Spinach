# examples/fundamentals/correlation_function_2.m

- MATLAB implementation: [examples/fundamentals/correlation_function_2.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/fundamentals/correlation_function_2.m)

- Signature: `correlation_function_2()`

## What is compared

The example compares a Monte Carlo estimate with Spinach's analytical expansion of `G(L,k,m,p,q)=<D{L}(k,m)*D{L}(p,q)'>`, using `k=-1, m=2, p=-1, q=2` at rank `L=2`. The model is axial rotational diffusion: the source sets `sigma_ax=0.1` and `sigma_eq=0.2`; the four Wigner indices are shifted by `L+1` for MATLAB array indexing.

## Numerical construction

With `1e6` steps and `nlags=300`, three standard-normal increments drive the direction-cosine matrix from the identity. In the source's explicit generator sum, `sigma_ax` weights the first skew-symmetric matrix (rotation about the z direction), while `sigma_eq` weights each of the other two. Every increment is applied by exponentiating the combined matrix and right-multiplying the current direction-cosine matrix. A `parfor` loop converts the trajectory to Euler angles and then to rank-`L` Wigner matrices. Normalised `xcorr` of the selected elements is scaled by `1/(2*L+1)`; after `ifftshift`, the first 300 entries pair with lag points 0 through 299.

## Analytical construction and observable

The dummy system uses zero magnet, isotope `G`, Redfield relaxation, lab-frame retention, zero equilibrium, spherical-tensor Liouville formalism, and no basis approximation. The source passes the two correlation times `1./(3*[sigma_ax sigma_eq].^2)` to the system. `corrfun(spin_system,L,k,m,p,q)` returns weights and rates used to build `sum_j weights{1}(j)*exp(rates{1}(j)*(0:(nlags-1)))`. The function plots the real Monte Carlo samples against that analytical curve; it defines no error tolerance or pass/fail assertion. The source estimates minutes of calculation time.
