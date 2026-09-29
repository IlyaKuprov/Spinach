# examples/esr_liq_pulsed/pulse_acquire_chrysene.m

[Source file](https://github.com/IlyaKuprov/Spinach/blob/main/examples/esr_liq_pulsed/pulse_acquire_chrysene.m)

## Call and result

Call **pulse_acquire_chrysene()** with no arguments. It returns no variables: the function builds local FID and spectrum variables and opens a figure using **kfigure** and **plot_1d**; no data file is written. The source labels calculation time as seconds.

## Spin system, symmetry, and relaxation

The example is described as W-band pulse-acquire FFT ESR of a chrysene cation radical in a non-viscous liquid. It imports **../standard_systems/chrysene_cation.log** with **gparse** and **g2spinach**, using label mapping **{{'E','E'},{'H','1H'}}** and passing **[0 0]**. **options.no_xyz=1** ignores coordinate information; the source comment says hyperfine couplings are provided. The magnetic-induction parameter is **3.5**.

Common-linewidth damping uses **'damp'**, diagonal retention, zero equilibrium, and damping value **1e6**. The **sphten-liouv** basis has no approximation, longitudinal **1H**, projection **+1**, and six **S2** symmetry pairs: **[1 7]**, **[2 8]**, **[3 9]**, **[4 10]**, **[5 11]**, and **[6 12]**.

## ESR acquisition

The initial state and receiver are **state(spin_system,'L+','E')**; **E** is detected and decoupling is empty. Offset is **-2e7**, sweep **1e8**, point count 1024, zero-fill 4096, axis label **'GHz-labframe'**, derivative 1, and axis inversion 1. The FID comes from **liquid(spin_system,@acquire,parameters,'esr')**; it receives **'none'** apodisation, is Fourier-transformed using the zero-fill length, and its real part is plotted.

## Dependencies and limits

Requires the relative standard-system log, Spinach **gparse**/**g2spinach** import helpers, and system/basis/state, liquid ESR/acquire, apodisation, FFT, and plotting routines. The source cites no paper DOI. It specifies no pulse shape or duration, so do not infer timing or RF details from the pulse-acquire label; the implemented acquisition is the **@acquire** call. The source values for field, damping, offset, and sweep carry no inline units; the axis-unit setting is explicit.
