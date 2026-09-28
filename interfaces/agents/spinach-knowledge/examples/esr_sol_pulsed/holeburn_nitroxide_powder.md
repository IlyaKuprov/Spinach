# examples/esr_sol_pulsed/holeburn_nitroxide_powder.m

- Signature: `holeburn_nitroxide_powder()`

## Purpose

Compares pulsed spectral hole burning for a nitroxide radical using a frequency-swept soft pulse, a fixed-frequency soft pulse, and a no-pulse reference. The soft-pulse simulations use the Fokker–Planck formalism; the script then acquires and Fourier transforms each signal.

## Model and calculation

The model contains an electron and `14N` at 3.5 T, with anisotropic electron Zeeman and electron–nitrogen coupling tensors. It uses a 64-point acquisition, zero-filled to 512 points, and a 3,200-point spherical powder grid. Each soft pulse has 100 one-nanosecond steps and rank 2: the swept pulse spans −350 to −250 MHz, while the fixed-frequency pulse is at −300 MHz. A third run sets the pulse power to zero as the reference. The signals are exponentially apodised and plotted as spectra. The source estimates a run time of seconds.
