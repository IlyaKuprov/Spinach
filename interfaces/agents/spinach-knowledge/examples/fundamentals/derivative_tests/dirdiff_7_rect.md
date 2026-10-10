# examples/fundamentals/derivative_tests/dirdiff_7_rect.m

- MATLAB implementation: [examples/fundamentals/derivative_tests/dirdiff_7_rect.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/fundamentals/derivative_tests/dirdiff_7_rect.m)

- Signature: `dirdiff_7_rect()`

## Question tested

For phase-modulated GRAPE with the rectangular integrator and configured waveform-distortion chains, does the analytical derivative of the first fidelity component agree with centred finite differences at selected phase samples?

## Setup

The script repeats the test for `sphten-liouv`, `zeeman-liouv`, and `zeeman-hilb`, using the system returned by `dirdiff_test_system`. Its controls use isotope `13C`, channel map `[1;1]`, drift `H`, operators `Lx` and `Ly`, initial states {Sx,Sy,Sz}, targets {-Sz,Sy,Sx}, power levels `2*pi*linspace(50e3,70e3,10)`, method `lbfgs`, maximum iterations 1000, and an empty plotting list. The rectangular grid and amplitude vector each have five entries: `12.8e-6*ones(1,5)` and `ones(1,5)`. The random phase guess is `randn(1,5)/3`; the increment is `h=sqrt(eps('double'))`.

The source configures two rows of four distortion callbacks:

1. `firf(w,[0.9 0.1i])`, `spf(w,0.2)`, `szf(w,0.2)`, then `amp_root(w,2*pi*20e3,4)`.
2. `szf(w,0.2)`, `spf(w,0.2)`, `amp_root(w,2*pi*20e3,4)`, then `firf(w,[0.9 0.1i])`.

These are the source-listed chains; no additional units are assigned here to their numerical arguments.

## Comparison and observable

The analytical phase gradient is returned by `grape_phase`. At phase entries 1, 3, and 5, the script evaluates the first fidelity entry for guesses perturbed by +h and -h, then forms `g_num=(fid_plus(1)-fid_minus(1))/(2*h)`. It compares this with the corresponding analytical-gradient entry using `abs(g_anl-g_num)/abs(g_num) < 1e-6`. Each result produces a formalism- and position-specific passed message or raises an error.

## Scope

The comparisons cover only the left edge, right edge, and middle entry of the five-element phase guess for each formalism. They do not establish agreement for all phase coordinates or report an outcome of a particular test run.