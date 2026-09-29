# examples/esr_liq_pulsed/pulse_acquire_methyl.m

[Source file](https://github.com/IlyaKuprov/Spinach/blob/main/examples/esr_liq_pulsed/pulse_acquire_methyl.m) · [Figure 4 reference](https://doi.org/10.1051/0004-6361:20020268)

## Call and result

Call **pulse_acquire_methyl()** with no arguments. It has no output arguments: it creates local FID and spectrum variables, then displays the real spectrum with **kfigure** and **plot_1d**; it does not write a data file. The source describes calculation time as seconds and identifies Figure 4 by Zhitnikov and Dmitriev as the target.

## Spin system and relaxation

The source describes an X-band pulse-acquire FFT ESR spectrum of methyl radical. It reads **../standard_systems/methyl.log** using **gparse** and **g2spinach**, maps **E** to **E** and **H** to **1H**, and passes **[0 0]**. With **options.no_xyz=1**, coordinate information is ignored; the source comment says hyperfine couplings are provided. It sets **sys.magnet=0.33**.

The simple common-linewidth model uses relaxation **'damp'**, diagonal retention, zero equilibrium, and **inter.damp_rate=2.5e7**. The **sphten-liouv** basis has no approximation, projection **+1**, and longitudinal **1H**; no symmetry group is specified.

## ESR acquisition

Detected spin is **E**; initial state and receiver are both **state(spin_system,'L+','E')**, with an empty decoupling list. Offset is **0**, sweep **5e8**, point count 256, zero-fill 1024, axis label **'GHz-labframe'**, derivative 1, and axis inversion 1. The function calls **liquid(spin_system,@acquire,parameters,'esr')**, applies **'none'** apodisation, Fourier-transforms using the zero-fill length, and plots the real spectrum.

## Dependencies and limits

Requires the relative **standard_systems/methyl.log** input, Spinach **gparse**/**g2spinach** import helpers, and system/basis/state, liquid ESR/acquire, apodisation, FFT, and plotting routines. The DOI cited by the source is linked above. The source values for field, damping, offset, and sweep have no units annotated in the file; the axis-unit field is explicitly **'GHz-labframe'**. No pulse shape or duration is specified: pulse-acquire is implemented through **@acquire**, not a pulse-program block in this function.
