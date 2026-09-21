# examples/esr_sol_pulsed/hpa_triplet.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/esr_sol_pulsed/hpa_triplet.m`
- Signature: `hpa_triplet()`
- Total lines: 73

## Purpose

Hypothetical powder averaged X-band pulse-acquire ESR spectrum of photogenerated pentacene triplet state. Calculation time: seconds.

## Physical / mathematical content

- Pulsed ESR / EPR solid-state examples. These scripts revolve around electron spin echo sequences, DEER, RIDME, ENDOR, ESEEM, and HYSCORE. They combine anisotropic Zeeman and hyperfine Hamiltonians with selective pulses, echo formation, and orientation averaging.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- Hypothetical powder averaged X-band pulse-acquire ESR
- spectrum of photogenerated pentacene triplet state.
- Calculation time: seconds.
- Magnet field
- Triplet electron
- Zeeman tensor, assumed isotropic
- ZFS, photo-excited pentacene triplet
- Basis set
- Disable trajectory-level SSR algorithms
- Spinach housekeeping
- Sequence parameters
- Zeeman tensor into Hz/Tesla

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `zfs2mat()`, `create()`, `basis()`, `state()`, `operator()`, `spin()`, `zftrip()`, `euler2dcm()`, `powder()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
