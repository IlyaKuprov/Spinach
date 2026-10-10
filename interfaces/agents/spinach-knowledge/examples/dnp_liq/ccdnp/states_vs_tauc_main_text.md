# examples/dnp_liq/ccdnp/states_vs_tauc_main_text.m

- MATLAB implementation: [examples/dnp_liq/ccdnp/states_vs_tauc_main_text.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/dnp_liq/ccdnp/states_vs_tauc_main_text.m)

- Signature: `states_vs_tauc_main_text()` (no arguments)
- Reference: [Journal of Magnetic Resonance 2021, 106940](https://doi.org/10.1016/j.jmr.2021.106940)
- The source comments that calculation takes seconds; this is not a benchmark guarantee.

## What it calculates

This main-text liquid-state DNP example follows the steady-state amplitudes of 12 selected spin observables across rotational correlation time for one `1H` coupled to two electrons. The electrons have anisotropic Zeeman tensors and scalar electron–electron exchange; the nucleus is coupled through the dipolar model represented by the supplied coordinates. A single field and microwave setting are held fixed during the sweep.

## Running and model inputs

Run in MATLAB with Spinach and its example/plotting functions on the MATLAB path. The function takes no arguments and loads no external data: the three-spin model and all sweep settings are defined in the source. It uses MATLAB `parfor`; the file does not configure a parallel pool.

- Field: 14.1 T. Isotope order is `{'1H','E','E'}` (nucleus, electron 1, electron 2).
- Nuclear Zeeman entries are `[0 10 20]` with Euler angles `[0 0 0]`. The electron tensors are `[2.0034 2.0038 2.0038]` at `[-0.872 -0.013 0.868]`, and `[2.0057 2.0030 2.0030]` at `[-1.145 0.061 1.143]`.
- Scalar coupling is set at `inter.coupling.scalar{2,3}=3e6`. Coordinates are `[0 0 0]`, `[5.090 0.010 0.958]`, and `[-5.090 0.061 1.032]`. The source does not state units for the coupling value, coordinates, or Euler angles.
- Use the full `sphten-liouv` basis (`bas.approximation='none'`) and Redfield relaxation with `inter.equilibrium='zero'`, secular retention, and `inter.temperature=298`. The source says the relaxation-integration tolerance `sys.tols.rlx_integration=1e-10` “needs to be this tight”; no temperature unit is stated there.
- The DNP settings are `parameters.spins={'E'}`, `mw_pwr=2*pi*1e6`, `mw_frq=2*pi*(-0.62e6)`, `method='lvn-backs'`, and `needs={'rho_eq'}`; `g_ref` is the mean of electron 1's listed Zeeman values. No unit is given for the microwave expressions in this source.

The 64-point correlation-time grid runs from `50e-12` to `500e-12` (the plot labels this axis in ps). At each point the script constructs the spin system and basis, builds 12 coil/detection states, and calls `liquid(spin_system,@dnp_freq_scan,locpar,'esr')` for the steady state.

## Output and limits

The three panels make the 12 plotted states explicit: panel 1 contains `E1+ + 2 E1+E2z`, `E1+ - 2 E1+E2z`, `E2+ + 2 E1zE2+`, and `E2+ - 2 E1zE2+`; panel 2 contains `2 NzE1+`, `2 NzE2+`, `4 NzE1+E2z`, and `4 NzE1zE2+`; panel 3 contains `2 NzE1z`, `2 NzE2z`, `4 NzE1zE2z`, and `Nz`. Each curve is the absolute value of its corresponding result entry, plotted as steady-state amplitude in arbitrary units against correlation time in ps. The function does not return the result array or save a data file; the source only creates the figure. It samples the specified single model and fixed microwave setting, not a frequency or field map. Numerical curves have not been reproduced here.
