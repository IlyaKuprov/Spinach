# examples/nmr_solids/mas_powder_suc_fplanck.m

Source: [examples/nmr_solids/mas_powder_suc_fplanck.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_solids/mas_powder_suc_fplanck.m)

## Purpose

The source header describes a 13C sucrose-powder MAS spectrum using the Fokker-Planck MAS formalism. It estimates “days” of calculation time; this is a source comment, not a measured timing. The filename and header identify the intended formalism, while the executable simulation call in this file is `singlerot`.

## Model and acquisition

The model is imported from `../standard_systems/sucrose.log` through `gparse` and `g2spinach`, selecting `13C` and using `182.1` as the absolute shielding reference (ppm; `g2spinach` defines the reference as placed at zero ppm). The source describes the chemical-shielding tensors, J-couplings, and coordinates as DFT-derived. The field is `14.1` T. The basis uses `sphten-liouv`, `IK-0`, projection `+1`, and inter-level 3; interaction and proximity cutoffs are 5.0 and 4.0.

The rotor-axis vector is `1 1 1`, MAS rate is `6000` Hz, maximum rank is 23, sweep is `50000` Hz, and offset is `15000` Hz. The acquisition has 256 points and zero-fills to 1024. It selects `13C`; `g2spinach` imports only carbon atoms, so proton spins and carbon–proton couplings are absent. This implements the header’s ideal `1H`-decoupling assumption without an RF channel; the script does not assign a `decouple` field. It assigns no dipolar tensor value directly; the parsed model also carries molecular coordinates. This source does not set `axis_units` or `invert_axis` explicitly.

## Calculation and display

Despite the filename and Fokker-Planck description, the code calls `singlerot(spin_system,@acquire,parameters,'nmr')`, with `parameters.grid='leb_2ang_rank_23'`; it does not call a function named `fplanck`. The initial state and receiver are both `L+` on `13C`. It exponentially apodises the calculated FID with parameter 6, Fourier transforms to 1024 points, and plots the real spectrum with `plot_1d`. The settings and resulting spectrum are a simulation, not an experimental measurement.
