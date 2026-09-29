# examples/esr_liq_pulsed/endor_nitroxide.m

- Call: `endor_nitroxide()` (no arguments; settings are defined in the function, with one external Gaussian-output dependency).
- Source: [`examples/esr_liq_pulsed/endor_nitroxide.m`](https://github.com/IlyaKuprov/Spinach/blob/main/examples/esr_liq_pulsed/endor_nitroxide.m)

## Spin system and input dependency

This is liquid-state Mims ENDOR for a `15N`-labelled nitroxide radical. Instead of defining its magnetic tensors inline, the function parses `../standard_systems/nitroxide.log` with `gparse`, then calls `g2spinach` using isotope mapping `{{'E','E'},{'N','15N'}}` and the `[0 0]` argument. It sets `options.no_xyz=1` because coordinates are not used, then sets `sys.magnet=0.33`. The field and imported magnetic values have no units stated in this source. Running the example therefore depends on that relative-path Gaussian output being available from the expected example working directory; the nitroxide parameters are not self-contained in this file.

## Simulation and output

The basis is `sphten-liouv` with `approximation='none'`; `sys.disable={'pt'}` disables path tracing. The function calls `liquid(spin_system,@endor_mims,parameters,'esr')` with offset 0, 512 points, sweep `1e8`, `tau=100e-9` (100 ns), zero filling to 4096, spin channel `{'E'}`, and axis units MHz. It mean-centres the FID, applies Kaiser apodisation with parameter 6, Fourier-transforms, and plots the magnitude with the x-axis labelled “Nuclear frequency, MHz”. It opens a figure but saves no data file. The source estimates a calculation time of seconds.

## Scope and caveats

The isotope mapping explicitly requests `15N`, while the numerical g and hyperfine parameters come from the external log rather than constants visible here. No explicit relaxation settings or DOI for the parameter calculation are supplied in the source.
