# examples/fundamentals/derivative_tests/dirdiff_4_trap.m

- MATLAB implementation: [examples/fundamentals/derivative_tests/dirdiff_4_trap.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/fundamentals/derivative_tests/dirdiff_4_trap.m)

- Signature: `dirdiff_4_trap()`

## Question tested

For phase-modulated GRAPE with the trapezium integrator, does the analytical derivative of the first fidelity component with respect to selected phase samples agree with a centred finite difference?

## Setup

The script repeats the comparison for `sphten-liouv`, `zeeman-liouv`, and `zeeman-hilb`, using the system returned by `dirdiff_test_system`. The control configuration uses isotope `13C`, channel map `[1;1]`, drift `H`, controls `Lx` and `Ly`, initial states {Sx,Sy,Sz}, and target states {-Sz,Sy,Sx}. It sets power levels to `2*pi*linspace(50e3,70e3,10)`, method `lbfgs`, maximum iterations 1000, and an empty plotting list. The trapezium grid has four entries, `12.8e-6*ones(1,4)`, while the amplitude vector and random phase guess each have five entries; the guess is `randn(1,5)/3`. The finite-difference increment is `h=sqrt(eps('double'))`.

## Comparison and observable

For phase entries 1, 3, and 5, the script forms `g_num=(fid_plus(1)-fid_minus(1))/(2*h)`, where `fid_plus` and `fid_minus` are the first fidelity entries returned by `grape_phase` at guesses differing by +h and -h in that entry. It compares this with the corresponding analytical-gradient entry from `grape_phase`. Each comparison uses the strict condition `abs(g_anl-g_num)/abs(g_num) < 1e-6`; the script prints a formalism-specific passed message or raises an error for that sample.

## Scope

Only the left edge, right edge, and middle sample of the five-element phase vector are checked in each formalism. This is not a check of every gradient coordinate or a reported result for a test run.