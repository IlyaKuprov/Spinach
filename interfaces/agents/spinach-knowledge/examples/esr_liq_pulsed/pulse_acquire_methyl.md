# examples/esr_liq_pulsed/pulse_acquire_methyl.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/esr_liq_pulsed/pulse_acquire_methyl.m`
- Signature: `pulse_acquire_methyl()`
- Total lines: 64

## Purpose

X-band pulse-acquire FFT ESR spectrum of methyl radical. Simple common line width is used as a relaxation model. Set to reprodu- ce Figure 4 from the paper by Zhitnikov and Dmitriev: Calculation time: seconds

## Physical / mathematical content

- Liquid-state ESR examples. The dominant physics is electron Zeeman interaction, hyperfine coupling, relaxation broadening, and pulse-acquire or ENDOR-type detection in fast tumbling systems.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- X-band pulse-acquire FFT ESR spectrum of methyl radical. Simple
- common line width is used as a relaxation model. Set to reprodu-
- ce Figure 4 from the paper by Zhitnikov and Dmitriev:
- Calculation time: seconds
- Ignore coordinate information (HFCs provided)
- Read the spin system (vacuum DFT calculation)
- Magnet induction
- Basis set
- Relaxation theory
- Spinach housekeeping
- Set the sequence parameters
- Simulation

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `g2spinach()`, `gparse()`, `create()`, `basis()`, `state()`, `liquid()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
