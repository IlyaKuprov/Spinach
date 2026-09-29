# examples/spin_chemistry/singlet_yield_4.m

## Purpose

The source comment associates this calculation with Figure 3 by Till, Timmel, Brocklehurst and Hore and cites [10.1016/S0009-2614(98)01158-0](https://doi.org/10.1016/S0009-2614(98)01158-0). The comment says the paper's field sweep uses only electron Zeeman operators and therefore omits the rising nuclear-Zeeman contribution on the high-field side. This is a report of what the source comment says, not an independent check of the cited paper. The source estimates calculation time as seconds; this is not a measured benchmark.

## Spin model and observable

The isotope list is two electrons and two protons (`{'E','E','1H','1H'}`). Their Zeeman scalar values are 2.0023, 2.0044, 0, and 0. The listed scalar hyperfine couplings connect electron 1 to the two protons at 35 G and 30 G, respectively. The source selects electron spins for the field-sweep Zeeman operator (`parameters.spins={'E'}` and `parameters.needs={'zeeman_op'}`), so the plotted sweep does not add nuclear Zeeman operators. The unapproximated basis is `zeeman-hilb`.

The shared `rydmr_exp` callback starts from the electron singlet selected by `[1 2]` and computes singlet recombination yield with exponential recombination kinetics. This script supplies no separate relaxation-rate parameter.

## Sweep and plot

The field array is `1e-3*10.^linspace(-5,3,2000)` T, spanning 10^-5 to 10^3 mT logarithmically. Five recombination rates are `[0.1 1.0 10.0 100.0 1000.0]*1e6` Hz. The plot uses `linspace(-5,3,2000)` as the log10 magnetic-induction axis in mT, and labels the vertical axis as singlet recombination yield.

Source: [examples/spin_chemistry/singlet_yield_4.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/spin_chemistry/singlet_yield_4.m).
