# examples/esr_liq_pulsed/pulse_acquire_benzoquinone.m

[Source file](https://github.com/IlyaKuprov/Spinach/blob/main/examples/esr_liq_pulsed/pulse_acquire_benzoquinone.m) · [Figure 1 reference](https://doi.org/10.1002/mrc.1260280313)

## Call and result

Call **pulse_acquire_benzoquinone()** with no arguments. It has no output arguments: the function builds local FID and spectrum variables, then opens a figure with **kfigure** and **plot_1d**. It does not save a data file. The source describes the calculation as taking seconds.

## Spin system and relaxation

The system is specified inline: one electron (**E**) and six protons (**1H**), with magnetic-induction setting **sys.magnet=0.33** and electron scalar Zeeman value **2.004577**. The six electron–proton scalar-coupling entries are **mt2hz(0.08)** three times, then **mt2hz(-0.059)**, **mt2hz(-0.364)**, and **mt2hz(-0.204)**. This represents the source-described 2-methoxy-1,4-benzoquinone radical in liquid state.

The common-linewidth relaxation model uses **'damp'**, diagonal retention, zero equilibrium, and **inter.damp_rate=1e6**. The basis is **sphten-liouv** with no approximation, longitudinal **1H**, projection **+1**, and **S3** symmetry over spins **[2 3 4]**.

## ESR acquisition

The initial state uses **state(spin_system,'L+','E')**; the detected spin is **E** and the decoupling list is empty. Offset is **-1e7**, sweep **3e7**, point count 1024, and zero-fill 4096. The axis label is **'GHz-labframe'**, and derivative and axis inversion are both 1. The function calls **liquid(spin_system,@acquire,parameters,'esr')**, applies **'none'** apodisation, computes **fftshift(fft(fid,parameters.zerofill))**, and plots the real spectrum. The receiver uses the same operator description with `coil_state` instead.

## Dependencies and limits

Requires Spinach system/basis/state construction, **mt2hz**, liquid ESR/acquire, apodisation, FFT, and plotting routines. Unlike the log-import examples, all spin-system parameters are supplied in this function; it reads no external spin-system file. It defines no explicit pulse shape, duration, or amplitude—the pulse-acquire calculation is represented by the **@acquire** liquid-ESR call. The field, damping, offset, sweep, and coupling values above are reproduced as coded; the source does not annotate their units. The axis-unit setting is explicit.
