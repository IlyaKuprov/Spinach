# examples/fundamentals/derivative_tests/dirdiff_6_rect.m

- MATLAB implementation: [examples/fundamentals/derivative_tests/dirdiff_6_rect.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/fundamentals/derivative_tests/dirdiff_6_rect.m)

- Signature: `dirdiff_6_rect()`

## Question tested

For Cartesian GRAPE with the rectangular integrator, do selected columns of the analytical Hessian agree with centred finite differences of the gradient?

## Setup

The script repeats the test for `sphten-liouv`, `zeeman-liouv`, and `zeeman-hilb`, using the system built by `dirdiff_test_system`. The control configuration uses isotope `13C`, channel map `[1;1]`, drift `H`, controls `Lx` and `Ly`, initial states {Sx,Sy,Sz}, targets {-Sz,Sy,Sx}, and power levels `2*pi*linspace(50e3,70e3,10)`. It sets a five-entry rectangular grid `12.8e-6*ones(1,5)`, method `newton`, `max_iter=1000`, an empty plotting list, guess `randn(2,5)/3`, and perturbation `h=sqrt(eps('double'))`.

## Comparison and observable

The analytical Hessian comes from `grape_xy`. For each selected linear index 1, 3, and 10, the script evaluates the Cartesian gradient at the guess perturbed by +h and -h at that index, flattens each gradient, and forms the numerical column `H_num=(grad_plus-grad_minus)/(2*h)`. It compares that column with the corresponding analytical Hessian column using the strict condition `norm(H_anl(:,j)-H_num,1) < 1e-6*norm(H_num,1)`. For each column, the script prints a formalism-specific passed message or raises an error.

## Scope and indexing

Only Hessian columns 1, 3, and 10 are checked. The source labels index 3 the middle column; it is MATLAB linear index 3 of the 2-by-5 waveform, not the central index of a ten-element vector. This selected-column check does not cover every entry of the Hessian.