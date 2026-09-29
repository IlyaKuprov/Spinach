# examples/relaxation_theory/aniso_diff_test_2.m

[Source](https://github.com/IlyaKuprov/Spinach/blob/main/examples/relaxation_theory/aniso_diff_test_2.m)

This example builds a Redfield relaxation superoperator for an anisotropically shielded 1H–13C pair undergoing anisotropic rotational diffusion. Its output is the full matrix R, displayed in MATLAB; it is not an experimental data fit or a simulated pulse-acquisition signal.

The field uses 2*pi*950.33e6/spin('1H'). The shielding principal values, in ppm, are [10 20 30] for proton and [40 50 60] for carbon, with Euler-angle rows [0 pi/4 0] for both. A scalar coupling of 145.0 Hz is specified. Unlike the companion geometry-based example, this source supplies no Cartesian coordinates, so it does not define a dipolar coupling. The rotational-diffusion eigenvalues are D=[2.16e8 2.35e8 7.45e8]; Spinach receives tau_c={1./(6*D)} as the associated correlation-time parameters. No spectral-density curve is calculated or plotted by the script.

The model requests Redfield relaxation, zero equilibrium and lab-frame retention, using sphten-liouv with no basis approximation. There is no explicit cross-correlation selection in the source, nor does it decompose the resulting matrix into named cross terms. The example creates the system, evaluates relaxation(spin_system), and prints full(R).
