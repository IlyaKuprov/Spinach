# examples/imaging/phase_encoding_2d.m

A simple phase-encoded 2D imaging example; the header estimates seconds of calculation time and credits Ahmed Allami and Ilya Kuprov. [Source](https://github.com/IlyaKuprov/Spinach/blob/main/examples/imaging/phase_encoding_2d.m)

## Model and sequence

The model is one `1H` spin with `sys.magnet=5.9` and zero scalar chemical shift. The source does not annotate a unit for the magnetic-induction value. Relaxation uses `t1_t2`, diagonal terms, zero equilibrium, and rate values `R1=30.0` and `R2=70.0` (no rate units are stated). Path tracing and Krylov propagation are disabled. The basis is `sphten-liouv` with no approximation.

The sequence calls `imaging(spin_system,@phase_enc_2d,parameters)`. It sets zero offset, `image_size=[101 105]`, readout-gradient amplitude `4.3e-3 T/m` and phase-encoding amplitude `3.8e-3 T/m`. Their duration values are `2e-3` and `1e-3`; the source does not attach units to the durations or to `t_echo=0.025`. The sample geometry is `[0.30 0.25]` with `[108 90]` points and `{'period',3}` differentiation.

Relaxation operators come from `rlx_t1_t2`; the spatial maps `R1Ph` and `R2Ph` are loaded from `../../etc/phantoms/letter_a.mat`. Uniform initial and coil phantoms use `Lz` and `L+` states.

## Output and interpretation

The returned image is displayed beside the two loaded relaxation maps using `mri_2d_plot`. The sequence's `image_size` and spatial `npts` are distinct configured arrays (`[101 105]` versus `[108 90]`); do not conflate image matrix size with the phantom grid. This example has no explicit post-call FFT or apodisation block.
