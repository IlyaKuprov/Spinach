# examples/spin_chemistry/singlet_yield_3.m

## Purpose

The source comment associates this calculation with Figure 1 by Timmel, Till, Brocklehurst, McLauchlan and Hore and cites [10.1080/00268979809483134](https://doi.org/10.1080/00268979809483134). This records the example's citation, not an independent verification of the paper or a measured yield. The source comment estimates calculation time as seconds; it is not a benchmark.

## Spin model and observable

The system is two electrons and one proton (`{'E','E','1H'}`). The electron Zeeman factors are 2.0023 and the proton factor is zero. A 20 G scalar hyperfine coupling connects electron 1 to the proton. In `rydmr_exp`, the selected electrons `[1 2]` define the singlet input; the callback calculates the singlet recombination yield with exponential recombination kinetics. No independent relaxation-rate setting appears in this script. The calculation uses the full, unapproximated `zeeman-hilb` basis.

## Sweep and plot

The script evaluates 200 fields from zero to `3 * 20 / 1e4` T (three times the 20 G hyperfine field), and nine rates given by `2*pi*[0.005 0.02 0.05 0.1 0.15 0.2 0.3 0.5 2.0]*gauss2mhz(20)*1e6`. The callback documents field in tesla and rates in hertz. The plot uses `omega/a` from 0 to 3 on the horizontal axis, singlet recombination yield vertically, and limits the vertical range to 0.2–1.0.

Source: [examples/spin_chemistry/singlet_yield_3.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/spin_chemistry/singlet_yield_3.m).
