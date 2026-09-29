# examples/optimal_control/case_studies/Skinner_JMR_2011/toroid_response_map.m

[Source](https://github.com/IlyaKuprov/Spinach/blob/main/examples/optimal_control/case_studies/Skinner_JMR_2011/toroid_response_map.m) · [Skinner et al., J. Magn. Reson. 209, 282–290 (2011)](https://doi.org/10.1016/j.jmr.2011.01.026)

Computes the paper's *rectangular-pulse benchmark* (Fig. 8A), not its optimised shaped pulse (Fig. 8B). The toroid has radii 1–6 mm and RF field proportional to inverse radius, so 25 kHz at the outer wall becomes 150 kHz next to the central conductor. For the shown call, a y-phase 3.74 µs hard pulse is evaluated over ±1.5 kHz offsets and 51 equally spaced radii; the detected transverse magnetisation is integrated with equal radial weighting rather than uniform weighting over RF amplitudes.

`[response,fig]=toroid_response_map([0;25e3],3.74e-6,linspace(-1.5e3,1.5e3,51),linspace(1e-3,6e-3,51))` returns the radius-by-offset x-magnetisation array and detected offset profile, and plots their three-dimensional surface and radial integral. All four inputs are required. This call numerically simulates the paper’s rectangular-pulse benchmark; it is not an independent closed-form check of Spinach or an experimental dataset.
