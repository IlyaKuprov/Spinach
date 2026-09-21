# examples/esr_sol_pulsed/hard_3_pulse_deer_gd_2.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/esr_sol_pulsed/hard_3_pulse_deer_gd_2.m`
- Signature: `hard_3_pulse_deer_gd_2()`
- Total lines: 91

## Purpose

Gadolinium(III) DEER experiment. The calculation is done by brute- force time propagation and powder averaging. Outermost ZFS transi- tion is excited by the probe pulse and the central transition is excited by the pump pulse. The pulses are assumed to be ideal. Note: gadolinium spin echo is very sharp and difficult to catch in simulations because they do not include zero-field splitting distributions found in experim

## Physical / mathematical content

- Pulsed ESR / EPR solid-state examples. These scripts revolve around electron spin echo sequences, DEER, RIDME, ENDOR, ESEEM, and HYSCORE. They combine anisotropic Zeeman and hyperfine Hamiltonians with selective pulses, echo formation, and orientation averaging.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- Gadolinium(III) DEER experiment. The calculation is done by brute-
- force time propagation and powder averaging. Outermost ZFS transi-
- tion is excited by the probe pulse and the central transition is
- excited by the pump pulse. The pulses are assumed to be ideal.
- Note: gadolinium spin echo is very sharp and difficult to catch in
- simulations because they do not include zero-field splitting
- distributions found in experimental systems.
- Calculation time: minutes.
- Spin system properties
- Basis set
- Spinach housekeeping
- Probe pulse operator -bottom transition

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `sind()`, `cosd()`, `create()`, `basis()`, `pauli()`, `speye()`, `state()`, `operator()`, `powder()`, `apodisation()`, `fftshift()`, `ft_axis()`, `kfigure()`, `scale_figure()`, `subplot()`, `kxlabel()`.
