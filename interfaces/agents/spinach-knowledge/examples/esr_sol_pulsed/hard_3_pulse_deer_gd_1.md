# examples/esr_sol_pulsed/hard_3_pulse_deer_gd_1.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/esr_sol_pulsed/hard_3_pulse_deer_gd_1.m`
- Signature: `hard_3_pulse_deer_gd_1()`
- Total lines: 94

## Purpose

Gadolinium(III) DEER experiment at W-band using ideal pulses. Set to reproduce Figure 2b from the paper by Otting and co-authors: The calculation is done by brute-force time propagation and grid pow- der averaging. Central transitions are used on both gadolinium ions. Note: gadolinium spin echo is very sharp and difficult to catch in simulations because they do not include zero-field splitting distributions found in 

## Physical / mathematical content

- Pulsed ESR / EPR solid-state examples. These scripts revolve around electron spin echo sequences, DEER, RIDME, ENDOR, ESEEM, and HYSCORE. They combine anisotropic Zeeman and hyperfine Hamiltonians with selective pulses, echo formation, and orientation averaging.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.
- The spin physics includes through-space magnetic dipole-dipole coupling, a rank-2 anisotropic interaction with strong orientation dependence and characteristic secular/non-secular structure.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- Gadolinium(III) DEER experiment at W-band using ideal pulses. Set to
- reproduce Figure 2b from the paper by Otting and co-authors:
- The calculation is done by brute-force time propagation and grid pow-
- der averaging. Central transitions are used on both gadolinium ions.
- Note: gadolinium spin echo is very sharp and difficult to catch in
- simulations because they do not include zero-field splitting
- distributions found in experimental systems.
- Note: flip-flop terms in the inter-electron dipolar interaction are
- switched off ('deer-zz') to mimic the effects of slightly dif-
- ferent pulse frequencies in the experiment.
- Calculation time: minutes.
- Spin system properties

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `operator()`, `powder()`, `apodisation()`, `fftshift()`, `ft_axis()`, `kfigure()`, `scale_figure()`, `subplot()`, `kxlabel()`, `ktitle()`.
