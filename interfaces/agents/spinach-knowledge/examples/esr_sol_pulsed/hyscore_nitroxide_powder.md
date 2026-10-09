# examples/esr_sol_pulsed/hyscore_nitroxide_powder.m

- Signature: `hyscore_nitroxide_powder()`
- Source: [`examples/esr_sol_pulsed/hyscore_nitroxide_powder.m`](https://github.com/IlyaKuprov/Spinach/blob/main/examples/esr_sol_pulsed/hyscore_nitroxide_powder.m)
- Sequence: [`experiments/esr_hyperfine/hyscore.m`](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/esr_hyperfine/hyscore.m)
- Reference: [Szosenfogel and Goldfarb (Figure 2a)](https://doi.org/10.1080/00268979809483260)

## Aim and spin interactions

The example is a powder-averaged, time-domain HYSCORE simulation for a `14N` nitroxide at 0.350 T; the source says it is set up for Figure 2a of the cited paper. The two spins are nitrogen-14 ((I=1)) and an electron (`E`). The electron has scalar (g=2.0000); the nitrogen Zeeman term is set to zero. Their isotropic hyperfine coupling is (5.0×10^6) Hz. The nitrogen quadrupole tensor is made by `eeqq2nqi(2.4e6,0.5,1,[0 0 0])`: the supplied quadrupole frequency is 2.4 MHz, asymmetry is 0.5, spin is 1, and its Euler angles are zero. The calculation uses the full sphten Liouville-space basis, without approximation; trajectory-level SSR and colorbar display are disabled.

## HYSCORE timing and output

The driver starts from electron `Lz`, detects electron `L+`, and requests a 20 MHz sweep, `tau = 136 ns`, 128 points per dimension, and the `rep_2ang_800pts_sph` powder grid. The `hyscore` sequence helper applies π/2–τ–π/2, selects zero-order electron coherence, evolves the indirect dimension, applies a π pulse, and forms the direct-dimension echo signal with the corresponding τ and detection-pulse operations. The 2D signal has its mean removed and cosine-apodisation is applied in both dimensions; it is zero-filled to 256×256 and transformed with a 2D FFT. The absolute spectrum is displayed as positive contours in MHz. The script plots but does not save the spectrum; its header estimates a calculation time of seconds.
