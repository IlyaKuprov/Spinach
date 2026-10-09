# examples/dnp_liq/ccdnp/tau_scan_main_text.m

- MATLAB implementation: [examples/dnp_liq/ccdnp/tau_scan_main_text.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/dnp_liq/ccdnp/tau_scan_main_text.m)

- Signature: `tau_scan_main_text()` (no arguments)
- Reference: [Journal of Magnetic Resonance 2021, 106940](https://doi.org/10.1016/j.jmr.2021.106940)
- The source comments that calculation takes seconds; this is not a benchmark guarantee.

## What it calculates

This main-text liquid-state DNP example makes a two-dimensional map of the steady-state proton response against microwave-frequency offset and rotational correlation time. The model is one `1H` plus two electrons with anisotropic Zeeman interactions, scalar electron–electron exchange, electron–nuclear dipolar couplings represented by the coordinates, and Redfield relaxation.

## Running and model inputs

Run in MATLAB with Spinach and its example/plotting functions on the MATLAB path. The no-argument function embeds the system and scan ranges and reads no external data. It uses `parfor`; the source does not configure a parallel pool.

- The field is 14.1 T; isotope order is `{'1H','E','E'}`. The nuclear Zeeman entries are `[0 10 20]` with Euler angles `[0 0 0]`.
- Electron tensors are `[2.0034 2.0038 2.0038]` at `[-0.872 -0.013 0.868]` and `[2.0057 2.0030 2.0030]` at `[-1.145 0.061 1.143]`. The pair's scalar coupling entry is `3e6`; coordinates are `[0 0 0]`, `[5.090 0.010 0.958]`, and `[-5.090 0.061 1.032]`. No units are stated for these coordinates, scalar entry, or Euler angles.
- Use `sphten-liouv` with `approximation='none'`; relaxation is Redfield with zero equilibrium polarisation, secular retention, temperature entry 298, and `sys.tols.rlx_integration=1e-10` (the source says this tolerance needs to be this tight). Temperature units are not specified in the file.
- DNP settings are `spins={'E'}`, `mw_pwr=2*pi*1e6`, `method='lvn-backs'`, `needs={'rho_eq'}`, with `g_ref` the mean of electron 1's Zeeman entries. Microwave offsets are set by `2*pi*linspace(-5,10,512)*1e6`; the figure labels this axis in MHz as offset from `g_iso^(1)`.

The scan uses 128 correlation times from `50e-12` to `500e-12`, plotted in ps. For each time, the script builds the spin system and basis, measures the proton `Lz` coil, sets electron microwave operators, and evaluates `liquid(spin_system,@dnp_freq_scan,localpar,'esr')` over the 512 offsets.

## Output and limits

The local answer is a complex 512-by-128 array (frequency by correlation time). At each correlation time it is divided by the isotropic-equilibrium expectation `localpar.coil'*rho_eq`, where `rho_eq=equilibrium(spin_system)`. The figure displays `real(answer)` as an image, with correlation time (ps) on x, microwave offset (MHz) on y, and a colorbar labelled steady-state `1H` DNP. The source does not save the map or return the array. It is one fixed field/system and a discrete grid; only the real normalised response is shown. Numerical map values have not been reproduced here.
