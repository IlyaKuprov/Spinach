# examples/esr_sol_pulsed/holeburn_gd_dota_powder.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/esr_sol_pulsed/holeburn_gd_dota_powder.m`
- Signature: `holeburn_gd_dota_powder()`
- Total lines: 93

## Purpose

A hole burning simulation for a gadolinium ion. The soft pulse is simulated using Fokker-Planck formalism. Zero-field splitting dis- ribution is sampled using the statistical parameters reported in Figure 5 of Raitsimring et al, App. Mag. Res. 28, 281-295 (2005). A numerical powder grid and numerical second-order rotating frame transformation are used. Note: non-central transition Gd(III) holes are very shallow. Calc

## Physical / mathematical content

- Pulsed ESR / EPR solid-state examples. These scripts revolve around electron spin echo sequences, DEER, RIDME, ENDOR, ESEEM, and HYSCORE. They combine anisotropic Zeeman and hyperfine Hamiltonians with selective pulses, echo formation, and orientation averaging.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- A hole burning simulation for a gadolinium ion. The soft pulse is
- simulated using Fokker-Planck formalism. Zero-field splitting dis-
- ribution is sampled using the statistical parameters reported in
- Figure 5 of Raitsimring et al, App. Mag. Res. 28, 281-295 (2005).
- A numerical powder grid and numerical second-order rotating frame
- transformation are used.
- Note: non-central transition Gd(III) holes are very shallow.
- Calculation time: minutes
- Initialize the spectra
- Get the sampling
- Get the figure going
- Loop over ZFS distribution

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `zfs_sampling()`, `kfigure()`, `zfs2mat()`, `create()`, `basis()`, `state()`, `powder()`, `apodisation()`, `fftshift()`, `plot_1d()`.
