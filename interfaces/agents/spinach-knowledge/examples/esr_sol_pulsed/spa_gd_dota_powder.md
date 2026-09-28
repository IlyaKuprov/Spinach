# examples/esr_sol_pulsed/spa_gd_dota_powder.m

- Signature: `spa_gd_dota_powder()`

## Purpose

Simulates a soft-pulse spectrum of a gadolinium ion using the Fokker–Planck formalism, powder averaging, and a third-order numerical rotating-frame transformation. The zero-field-splitting (ZFS) distribution is sampled using statistical parameters reported in Figure 5 of Raitsimring et al., App. Mag. Res. 28, 281–295 (2005). Calculation time: hours.

## Physical / mathematical content

- Models an `E8` electron spin with `sys.magnet=3.5`, scalar Zeeman parameter `2.002319`, and a ZFS tensor constructed from each sampled pair of `D` and `E` values.
- Uses a spherical powder grid (`rep_2ang_400pts_sph`), an `Lz` initial state, and an `L+` detection state. The rotating-frame setting is `{{'E8',3}}`.
- Applies a rank-2 soft pulse with phase `-pi/2`, frequency `-0.5e9`, duration `50 ns`, and power `2*pi*0.02e9`.

## Numerical / algorithmic content

- Obtains ZFS samples and weights with `zfs_sampling(30,5,1e-4)` and runs `powder(spin_system,@sp_acquire,parameters,'labframe')` for each sample.
- Applies exponential apodisation to each acquired FID, Fourier-transforms it with 2048-point zero filling, and adds the result to the spectrum weighted by its ZFS sampling weight.
- Uses a sweep of `0.8e10`, 512 acquisition points, a GHz plot axis, and the `expm` propagation method.

## Implementation structure

1. Preallocates a complex 2048-point spectrum, obtains the ZFS samples, and opens a figure.
2. For each sample, builds a spherical-tensor Liouville-space spin system without basis approximation, sets the acquisition and soft-pulse parameters, and simulates powder-averaged acquisition.
3. Apodises and Fourier-transforms the FID, accumulates the weighted spectrum, and plots its real part after each iteration.