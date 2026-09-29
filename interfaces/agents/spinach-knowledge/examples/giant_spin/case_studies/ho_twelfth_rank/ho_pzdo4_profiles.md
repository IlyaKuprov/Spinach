# examples/giant_spin/case_studies/ho_twelfth_rank/ho_pzdo4_profiles.m

- MATLAB implementation: [examples/giant_spin/case_studies/ho_twelfth_rank/ho_pzdo4_profiles.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/giant_spin/case_studies/ho_twelfth_rank/ho_pzdo4_profiles.m)

- Signature: `ho_pzdo4_profiles()`

## Purpose

Calculates pulsed-field magnetisation of the Ho(pzdo)4 metal-organic framework, represented by a J=8 giant spin with crystal-field terms through spherical rank 12, for four field profiles. The case study corresponds to Fig. 2 of [arXiv:2609.16352](https://arxiv.org/abs/2609.16352).

## Physical model and setup

The system uses one `E17` effective spin with g = 1.24 and the CASSCF crystal-field coefficients from `ho_pzdo4_params.m`, converted rank by rank from cm^-1 to Hz and then to spherical form. The single-crystal field frame is aligned with the laboratory frame (orientation [0 0 0]). The formalism is Zeeman-Hilbert space without approximation. Spin-phonon relaxation uses the generalised Lindblad dissipator of Saito and Miyashita with a super-Ohmic bath (alpha = 2, lambda = 10 cm^-1, I0 = 1e-14 ps/rad); the bath parameter is converted to rad/s units in the source. The coupling matrix has unit elements between adjacent mJ states. Magnetisation is plotted as the moment -1.24 Jz in Bohr magnetons.

## Field profiles and propagation

Field is in tesla and time in seconds in the profile functions:

1. Linear sweep B(t) = 1e4 t, i.e. 10 T/ms.
2. Piecewise-linear interpolation through (0, 0), (1 us, 0.1 T), (10 us, 1 T), (100 us, 5 T), and (1 ms, 10 T).
3. Monotone cubic interpolation (`pchip`) through 76 measured points of the 65 T pulse, tabulated from 0 to 10 ms. The plotted propagation covers its first millisecond.
4. A 0.1 T sinusoid at angular frequency 0.134124264765e12 rad/s (about 0.71 cm^-1), corresponding to a period of about 46.8 ps.

The first three profiles are propagated at 2 K with 100,000 steps of 10 ns and output every 100 steps (`nout=100`), covering 1 ms. The sinusoidal profile uses 1 mK, 4,680 steps of 0.1 ps, and output every four steps (`nout=4`), covering 468 ps (ten periods).

## Output and plot limits

The function produces a four-panel plot: magnetisation is red on the left axis and field is black on the right. The first three panels use time limits 0-1e9 ps and magnetisation limits 0-6 Bohr magnetons; the sinusoidal panel uses 0-468 ps and -6 to 6 Bohr magnetons. Field limits are 0-10 T for the first two profiles, 0-13 T for the measured-pulse segment, and -0.25 to 0.25 T for the sinusoid. The plotted temperature and time step are annotated in each panel. These are the source's configured intervals and axes.
