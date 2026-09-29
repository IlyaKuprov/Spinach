# examples/esr_liq_pulsed/pulse_acquire_biaryl.m

[Source file](https://github.com/IlyaKuprov/Spinach/blob/main/examples/esr_liq_pulsed/pulse_acquire_biaryl.m)

## Call and result

Call **pulse_acquire_biaryl()** with no arguments. There are no returned values: the function constructs a local FID and spectrum, then displays the real spectrum through **kfigure** and **plot_1d**; it does not save a data file. The source estimates calculation time as seconds.

## Spin system and relaxation

The source describes a time-domain pulse-acquire version of the EasySpin biaryl test file and acknowledges Stefan Stoll. It defines one electron, two **14N** nuclei, and ten **1H** nuclei, with **sys.magnet=0.33** and electron scalar Zeeman value **2.00316**. Each of the two equivalent nitrogen/proton groups has the same electron-coupling entries: **12.16e6** for nitrogen and, for the five proton entries, **-6.7e6**, **-1.82e6**, **-7.88e6**, **-0.64e6**, and **67.93e6**. These are the raw values assigned in **inter.coupling.scalar**; the source does not annotate units.

Relaxation is **'damp'** with diagonal retention, zero equilibrium, and damping value **5e5**. The **sphten-liouv** basis uses no approximation, longitudinal **1H** and **14N**, projection **+1**, and six **S2** pairs: **[2 8]**, **[3 9]**, **[4 10]**, **[5 11]**, **[6 12]**, and **[7 13]**.

## ESR acquisition

The initial state and receiver are both **state(spin_system,'L+','E')**; detected spin is **E** and decoupling is empty. Acquisition settings are offset **0**, sweep **3e8**, 4096 points, zero-fill 16384, axis label **'GHz-labframe'**, derivative 1, and axis inversion 1. The function calls **liquid(spin_system,@acquire,parameters,'esr')**, applies no apodisation, Fourier-transforms the FID with the configured zero-fill, and plots its real part.

## Dependencies and limits

Requires Spinach system/basis/state, liquid ESR/acquire, apodisation, FFT, and plotting routines; all spin-system values are in the function, with no external log input. The source describes explicit time propagation in Liouville space and full **S2xS2xS2xS2xS2xS2** direct-product symmetry. It supplies no DOI or direct URL for the EasySpin test it references. No pulse shape or duration is defined; the implemented sequence path is **@acquire** through liquid ESR. Numeric settings are transcribed from the source without inferred physical units.
