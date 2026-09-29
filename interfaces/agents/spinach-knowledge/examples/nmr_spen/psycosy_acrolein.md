# examples/nmr_spen/psycosy_acrolein.m

- MATLAB implementation: [examples/nmr_spen/psycosy_acrolein.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_spen/psycosy_acrolein.m)

- Signature: `psycosy_acrolein()`
- Source: [`examples/nmr_spen/psycosy_acrolein.m`](../../../../../examples/nmr_spen/psycosy_acrolein.m)

## Spin model and spatial encoding

The model contains five protons at 14.1 T. The source assigns shift entries [9.11, 5.91, 5.49, 5.32, 7.32] and pairwise scalar-coupling entries (1,2)=7.7, (2,3)=10.0, (2,4)=17.3, (1,3)=1.0, and (1,4)=1.0; it also sets two entries involving site 5 to zero. The entries are stated as supplied by the source, without adding units not given there. The basis uses the spherical-tensor Liouville formalism, IK-2 approximation, proximal level 1, and scalar-coupling connectivity.

The `@psycosy` sequence is simulated with `imaging` on a 15 mm sample with 100 spatial points and derivative option `{'period',3}`. Diffusion and flow are zero. A zero-valued relaxation phantom and the relaxation superoperator are supplied, with `Lz` initial and `L+` detection states.

## Sequence and processing

The sequence parameters are offset 4392, sweep 3000, [512, 512] acquisition points, [1024, 1024] zero-fill points, and `axis_units='ppm'`. The mixing time is 70e-3 s and gradient amplitude 1e-2 T/m. The saltire chirp uses a 20-degree flip angle, 0.015 s pulse duration, 0.05 s chirp-gradient duration, and 10000 Hz sweep width. The source does not annotate units for the offset or the `sweep` assignment.

The simulated 2D FID is square-sine apodised in both dimensions and transformed with a zero-filled 2D FFT. The plot displays the magnitude of the spectrum. The source header estimates hours of calculation and notes faster execution on a GPU; GPU enablement is commented out.
