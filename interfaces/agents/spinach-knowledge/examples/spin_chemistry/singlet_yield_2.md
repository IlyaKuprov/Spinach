# examples/spin_chemistry/singlet_yield_2.m

## Purpose

This example scans magnetic field and singlet recombination rate for a liquid-state radical pair with six equivalent protons. Its source comment describes an exponential recombination-kinetics model and full S6 symmetry, and estimates the example calculation time as seconds; that timing is a source comment, not a measured benchmark.

## Spin model and observable

The isotope list is two electrons and six protons (`{'E','E','1H','1H','1H','1H','1H','1H'}`). Both electron Zeeman factors are 2.002; nuclear Zeeman factors are set to zero. The coupling matrix connects electron 1 to each of the six protons with the same 0.295 mT scalar hyperfine value; electron 2 has no listed hyperfine coupling. The basis uses `sphten-liouv`, no approximation, projection 0, and full S6 permutation symmetry over spin indices 3–8.

The shared `rydmr_exp` callback forms an electron singlet on the selected pair `[1 2]` and calculates singlet recombination yield. The script supplies no separate relaxation-rate parameter.

## Sweep and plot

`parameters.fields` runs from 0 to 5 mT in 0.01 mT steps (501 field values, passed to the callback in tesla). The recombination-rate list is `[0.176 0.880 1.76 3.52 8.8 17.6 35.2 52.8] * 1e6` Hz. The plotted matrix `M` is the singlet yield versus field for these rate values; the field axis is labelled in tesla.

Source: [examples/spin_chemistry/singlet_yield_2.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/spin_chemistry/singlet_yield_2.m).
