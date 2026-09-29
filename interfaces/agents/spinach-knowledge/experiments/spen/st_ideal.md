# experiments/spen/st_ideal.m

[Canonical source](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/spen/st_ideal.m) · [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=st_ideal.m)

Ideal Stejskal–Tanner sequence in the notation of Figure 1 in [the cited article](http://dx.doi.org/10.1002/cmr.a.21241), with no gaps between sequence events. The source forms L = H + F + 1i*R + 1i*K and applies a pi/2 pulse about Ly to rho0. It evolves under L+g_amp*G{1} for delta_sml, evolves freely for (delta_big-delta_sml)/2, applies a pi pulse about Ly, repeats the free delay, and then applies the second L+g_amp*G{1} interval for delta_sml. The output inten = abs(coil'*rho) is a scalar: the absolute detected first FID point, not a sampled multidimensional FID. The source describes this value as proportional to the integral of the real part of the correctly phased spectrum.

Call from the imaging() context, which supplies H, R, K, G, and F. Required parameters fields are rho0, coil, scalar npts (spin-packet count), g_amp (T/m), positive finite real scalar delta_sml, positive finite real scalar delta_big no smaller than delta_sml, and a single working-spin entry in spins. The function requires sphten-liouv; H, R, K, and F must be equal-sized matrices and G must be a cell array. The timing quantities are passed directly to the source propagators; the source header gives the gradient-amplitude unit but does not state an explicit time unit.

The source credits mariagrazia.concilio@sjtu.edu.cn and ilya.kuprov@weizmann.ac.il.
