# examples/esr_liq_pulsed/endor_methyl.m

- Call: `endor_methyl()` (no arguments; settings are defined in the function).
- Source: [`examples/esr_liq_pulsed/endor_methyl.m`](https://github.com/IlyaKuprov/Spinach/blob/main/examples/esr_liq_pulsed/endor_methyl.m)

## Spin system

The example simulates liquid-state Mims ENDOR of a methyl radical. Its source describes the magnetic parameters as DFT-derived, but does not identify a calculation or external input file. The model contains three equivalent protons and one electron: `sys.isotopes={'1H','1H','1H','E'}`, electron Zeeman scalar 2.0026, and three equal electron–proton scalar entries of `-64.4e6`. The source does not state units for the field value `sys.magnet=0.346` or those raw coupling values. The three proton spins are treated with `S3` symmetry; the full `sphten-liouv` basis is configured with `approximation='none'`.

## Simulation and output

The call is `liquid(spin_system,@endor_mims,parameters,'esr')`. Parameters set offset 0, 512 points, sweep `1.2e8`, `tau=100e-9` (100 ns), zero filling to 4096, detected spin `{'E'}`, and axis units MHz. The example subtracts the FID mean, applies Kaiser apodisation with parameter 6, Fourier-transforms, and plots the spectrum magnitude. Its x-axis label is explicitly “Nuclear frequency, MHz”. It creates a figure but writes no data file; the source estimates a calculation time of seconds.

## Scope and caveats

The `spins={'E'}` setting is the channel supplied to the Mims ENDOR simulation, while the plotted label says nuclear frequency; preserve both source settings rather than interpreting one as a typographical error. The example contains no explicit relaxation configuration and gives no DOI for the stated DFT-derived parameters.
