# examples/esr_sol_pulsed/spa_gd_dota_powder.m

- Function: `spa_gd_dota_powder()`.

## Model

This example builds a powder-averaged soft-pulse EPR calculation for a gadolinium ion with sampled zero-field splitting (ZFS). Its source comment attributes the distribution's statistical parameters to Figure 5 of Raitsimring et al., *App. Mag. Res.* 28, 281-295 (2005), and describes a third-order numerical rotating-frame transformation. The source estimates hours. It calls `[D,E,W]=zfs_sampling(30,5,1e-4)` and loops over the returned samples and weights; the file does not state units for these sampling arguments.

For each sample the code sets `sys.magnet=3.5`, isotope `E8`, Zeeman scalar `2.002319`, and self-coupling matrix `0.56e9*zfs2mat(D(n),E(n),0,0,0)`. The basis is `sphten-liouv` with no approximation and projections `-3:3`; trajectory-level SSR is disabled. The rotating-frame setting is `parameters.rframes={{'E8',3}}`. The source does not annotate units for its field or coupling values.

## Pulse, acquisition, and plotted observable

The soft pulse has rank 2, phase `-pi/2`, frequency `-0.5e9`, duration `50.0e-9` s (50 ns), power `2*pi*0.02e+9`, and method `expm`. Raw frequency and power units are not annotated. Initial state and receiver are `Lz` and `L+` on `E8`; no spins are listed for decoupling. Acquisition settings are offset 0, sweep `0.8e10`, 512 points, zero-fill 2048, `axis_units='GHz'`, grid `rep_2ang_400pts_sph`, derivative off, and axis inversion off.

For each sample, `powder(spin_system,@sp_acquire,parameters,'labframe')` acquires an FID. The script applies exponential apodisation with parameter 10, computes a 2048-point shifted FFT, and adds `W(n)` times that transform to a complex spectrum accumulator. It plots `real(spectrum)` with `plot_1d` inside the sample loop, so the displayed sum is progressively accumulated. There is no explicit final normalisation in this script.

## Source

[examples/esr_sol_pulsed/spa_gd_dota_powder.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/esr_sol_pulsed/spa_gd_dota_powder.m)
