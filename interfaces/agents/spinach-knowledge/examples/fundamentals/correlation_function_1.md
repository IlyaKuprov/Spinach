# examples/fundamentals/correlation_function_1.m

- MATLAB implementation: [examples/fundamentals/correlation_function_1.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/fundamentals/correlation_function_1.m)

- Signature: `correlation_function_1()`

## What is compared

This example compares a finite Monte Carlo trajectory with Spinach's analytical correlation-function expansion for `G(L,k,m,p,q)=<D{L}(k,m)*D{L}(p,q)'>`. The source uses the prime on the second Wigner element and selects `k=-1, m=2, p=-1, q=2`; it converts each index from `[-L,L]` to MATLAB indexing by adding `L+1`. The case is isotropic rotational diffusion, with `sigma_iso=0.1` and rank `L=2`.

## Numerical construction

The Monte Carlo trajectory has `1e6` steps and requests `nlags=300`. At each step, three independent standard-normal increments weight the three displayed skew-symmetric rotation generators equally by `sigma_iso`; their sum is exponentiated and right-multiplied into the direction-cosine matrix, starting from `eye(3)`. A `parfor` loop converts each saved matrix through `dcm2euler` and `wigner(L,...)`; the selected Wigner elements are passed to normalised `xcorr`. The code multiplies that result by `1/(2*L+1)`, applies `ifftshift` to the correlation and lag vector, and plots the first 300 entries (lag points 0 through 299).

## Analytical construction and observable

A one-spin dummy system (zero magnet, isotope `G`) is configured with Redfield relaxation, lab-frame retention, zero equilibrium, spherical-tensor Liouville formalism, and no basis approximation. Its correlation time is `1/(3*sigma_iso^2)`. After `create` and `basis`, `corrfun(spin_system,L,k,m,p,q)` supplies weights and rates; the curve is assembled as `sum_j weights{1}(j)*exp(rates{1}(j)*(0:(nlags-1)))`. The function displays real Monte Carlo points against this curve. It has no numerical tolerance or pass/fail assertion; the plotted comparison, not an automated verdict, is the observable. The source estimates minutes of calculation time.
