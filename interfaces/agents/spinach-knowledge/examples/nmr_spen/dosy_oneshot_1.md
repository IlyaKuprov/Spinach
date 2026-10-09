# examples/nmr_spen/dosy_oneshot_1.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_spen/dosy_oneshot_1.m)

- Signature: dosy_oneshot_1()

## Purpose

Simulates a one-shot DOSY imaging sequence for three coupled 1H spins with different longitudinal and transverse relaxation rates. The source estimates minutes on an NVIDIA Tesla A100 and substantially longer on CPU.

## Spin system, relaxation, and acquisition

At a field parameter of 11.7428 T, the three shift values are 4.70, 3.50, and 1.50. The scalar-coupling network has pair values 15 for spins 1–2, 15 for 2–3, and 10 for 1–3. The t1_t2 relaxation model keeps diagonal relaxation terms and uses zero equilibrium; the reciprocal-rate inputs are formed from T1 values [0.1952, 0.2100, 0.2500] and T2 values [0.1602, 0.1802, 0.1902].

The acquisition sweep is 5000 Hz with 1024 points, zero-filled to 32768; the axis is ppm and the offset is 2497.78 Hz. The 0.015 m sample is represented by 5000 spatial points with a 7-point periodic derivative stencil. Uniform spatial phantoms set the initial Lz state and L- detection state.

## Diffusion encoding and observable

The reference diffusion coefficient is 18.55 × 10⁻¹⁰ m²/s. The dosy_oneshot imaging sequence uses gradient amplitude 0.255 T/m, kappa 0.2, gradient duration 0.001 seconds, diffusion delay 0.05 seconds, and gradient-stabilisation delay 0.0005 seconds. The resulting FID is exponentially apodised with parameter 5, Fourier transformed with the specified zero filling, and plotted as the negative real spectrum.
