# examples/spin_chemistry/cidnp_pumping_2.m

- Signature: `cidnp_pumping_2()`
- Source: [examples/spin_chemistry/cidnp_pumping_2.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/spin_chemistry/cidnp_pumping_2.m)

## Purpose

Evolves three coupled proton/fluorine detection channels with relaxation and pumping terms. The source header identifies this as a simulation of Figure 2A in the author's chemically amplified-NOE paper and points to [DOI 10.1016/j.jmr.2004.01.011](https://doi.org/10.1016/j.jmr.2004.01.011); the DOI is a source-comment citation, not independent confirmation of the plotted result. The script does not define radical-pair recombination kinetics.

## Spin system and evolution

- The two spins are `1H` and `19F` at `sys.magnet=14.1` T. The shifts are zero; the source labels the fluorine CSA values `[-47 -16 63]` as DFT, the coordinates `[0 0 0]` and `[0 2.60 0]` as DFT, and the scalar coupling `50` Hz as experimental. The source does not specify units for the CSA entries or coordinates.
- It builds a Redfield relaxation matrix with secular terms, `tau_c=110e-12` s, equilibrium `IME` and temperature value `298`; the Hamiltonian is obtained with the `nmr` assumption. The initial operator is `unit_state + 2*Hz + 2*Fz`.
- Pumping terms are applied to `2*Hz` and `2*Fz` with strengths `1.3` and `34.0`, respectively, and an additional proton relaxation term `-3.0*(Hz*Hz')` is subtracted from `R`. The three detected channels are `[Fz -HzFz Hz]`, propagated with `evolution(...,0.1,40,'multichannel')`.

## Visualisation

The figure plots the real parts of the three output channels, labelled `Fz`, `-HzFz` and `Hz`, against seconds. The plot's `time_axis` is `linspace(0,4,41)`. In `evolution`, the corresponding inputs are a timestep of `0.1` s and `40` steps, so the 41 plotted samples span 0–4 s; this agrees with the API signature in `kernel/evolution.m`. The function reports calculation time as seconds.
