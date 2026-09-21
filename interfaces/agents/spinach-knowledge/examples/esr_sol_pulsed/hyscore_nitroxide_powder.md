# examples/esr_sol_pulsed/hyscore_nitroxide_powder.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/esr_sol_pulsed/hyscore_nitroxide_powder.m`
- Signature: `hyscore_nitroxide_powder()`
- Total lines: 65

## Purpose

Powder-averaged HYSCORE on a 14N nitroxide radical. Time-domain simulation in Liouville space. Set to reproduce Figure 2a from Calculation time: seconds

## Physical / mathematical content

- Pulsed ESR / EPR solid-state examples. These scripts revolve around electron spin echo sequences, DEER, RIDME, ENDOR, ESEEM, and HYSCORE. They combine anisotropic Zeeman and hyperfine Hamiltonians with selective pulses, echo formation, and orientation averaging.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- Powder-averaged HYSCORE on a 14N nitroxide radical. Time-domain
- simulation in Liouville space. Set to reproduce Figure 2a from
- Calculation time: seconds
- Magnet field
- System specification
- Basis set
- Disable trajectory-level SSR algorithms
- Spinach housekeeping
- Set the sequence parameters
- Simulation
- Centre signal suppression
- Apodisation

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `eeqq2nqi()`, `create()`, `basis()`, `state()`, `powder()`, `apodisation()`, `fftshift()`, `fft2()`, `kfigure()`, `scale_figure()`, `plot_2d()`.
