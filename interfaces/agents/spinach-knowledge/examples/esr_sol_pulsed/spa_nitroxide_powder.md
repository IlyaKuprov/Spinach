# examples/esr_sol_pulsed/spa_nitroxide_powder.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/esr_sol_pulsed/spa_nitroxide_powder.m`
- Signature: `spa_nitroxide_powder()`
- Total lines: 73

## Purpose

A soft pulse simulation for a nitroxide radical powder. The soft pulse is simulated using the Fokker-Planck formalism; it is fol- lowed by time domain acquisition and Fourier transform. Calculation time: seconds

## Physical / mathematical content

- Pulsed ESR / EPR solid-state examples. These scripts revolve around electron spin echo sequences, DEER, RIDME, ENDOR, ESEEM, and HYSCORE. They combine anisotropic Zeeman and hyperfine Hamiltonians with selective pulses, echo formation, and orientation averaging.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- A soft pulse simulation for a nitroxide radical powder. The soft
- pulse is simulated using the Fokker-Planck formalism; it is fol-
- lowed by time domain acquisition and Fourier transform.
- Calculation time: seconds
- Isotopes
- Magnet field
- Interactions
- Basis set
- Disable trajectory-level SSR algorithms
- Spinach housekeeping
- Sequence parameters
- Soft pulse parameters

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `powder()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
