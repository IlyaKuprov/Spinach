# examples/esr_sol_pulsed/spa_nitroxide_powder.m

- Signature: `spa_nitroxide_powder()`

## Purpose

Simulate a soft-pulse spectrum of a nitroxide radical powder using the Fokker–Planck formalism, followed by time-domain acquisition and Fourier transformation. The source notes a calculation time of seconds.

## Physical / mathematical content

- The spin system contains an electron and `14N` at a magnetic field of 3.5 T. It specifies an anisotropic electron Zeeman matrix and an electron–nitrogen coupling matrix, including off-diagonal x–z terms.
- The calculation uses the `sphten-liouv` basis without approximation and the `rep_2ang_3200pts_sph` powder grid.
- The initial state is electron `Lz`; the detection state is electron `L+`. No spins are decoupled.

## Numerical / algorithmic content

- The soft pulse has rank 2, phase `-pi/2`, frequency `-300e6` Hz, duration `100e-9` s, and power `2*pi*16.5e6`. Its propagation method is `expm`.
- Acquisition uses an offset of `-2e8` Hz, a sweep of `8e8` Hz, and 64 points. The FID receives `crisp` apodisation before an FFT zero-filled to 512 points and shifted with `fftshift`.

## Implementation structure

- Construct the spin system with `create`, set its basis with `basis`, and disable trajectory-level SSR algorithms through `sys.disable={'trajlevel'}`.
- Run `powder(spin_system,@sp_acquire,parameters,'esr')` to acquire the powder-averaged FID.
- Plot the real spectrum with `plot_1d`, using MHz axis units.