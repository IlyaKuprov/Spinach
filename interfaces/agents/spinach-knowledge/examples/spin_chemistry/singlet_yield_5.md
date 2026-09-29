# examples/spin_chemistry/singlet_yield_5.m

## Purpose

The source comment associates this calculation with Figure 5 by Timmel, Till, Brocklehurst, McLauchlan and Hore and cites [10.1080/00268979809483134](https://doi.org/10.1080/00268979809483134). This preserves the script's citation; it does not independently verify the paper or report a measured yield. The source comment estimates calculation time as seconds, not as a benchmark.

## Spin model and observable

The system contains two electrons and two protons (`{'E','E','1H','1H'}`). Both electron Zeeman factors are 2.0023 and the nuclear factors are zero. The scalar hyperfine couplings are 20 G from electron 1 to proton 1 and 20 G from electron 2 to proton 2. The source uses an unapproximated `zeeman-hilb` basis. Its `rydmr_exp` callback starts from the singlet of electrons `[1 2]` and returns singlet recombination yield using exponential recombination kinetics. The script declares no separate relaxation-rate parameter.

## Sweep and plot

Two hundred fields cover zero to `0.5 * 20 / 1e4` T (1 mT). Seven rates are defined by `2*pi*[1e-4 1e-3 0.01 0.02 0.05 0.1 0.2]*gauss2mhz(20)*1e6`; the callback documents rate values in hertz. The plot uses `omega/a` from 0 to 0.5 horizontally, shows the singlet-yield output, and clips the vertical axis to 0.46–0.58.

Source: [examples/spin_chemistry/singlet_yield_5.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/spin_chemistry/singlet_yield_5.m).
