# examples/relaxation_theory/dd_csa_xcorr_1.m

- Signature: `dd_csa_xcorr_1()`
- Source: [examples/relaxation_theory/dd_csa_xcorr_1.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/relaxation_theory/dd_csa_xcorr_1.m)

This example builds a BRW relaxation superoperator for anisotropically shielded `1H` and `13C` spins with a through-space dipole–dipole interaction. The source identifies both CSA–CSA and dipole–CSA (DD–CSA) cross-correlations and says the dipolar coupling is calculated from Cartesian coordinates. Spinach forms the relaxation contributions from the specified interactions; the example does not select individual cross terms by hand.

At `14.1 T`, the shielding principal values are `[7 15 -22]` and `[11 18 -29]` ppm, oriented by Euler triplets `[pi/3 pi/4 pi/5]` and `[pi/6 pi/7 pi/8]` radians. The coordinates are `[0 0 0]` and `[0 0 1.02]` Å, a separation of `1.02 Å` along the coordinate z direction. They supply the geometry for the dipolar interaction. Redfield relaxation uses one correlation time, `tau_c=1e-9 s`, zero equilibrium, and `rlx_keep='labframe'`; the basis is `sphten-liouv` with no approximation. A single correlation time describes one rotational-correlation timescale rather than an explicitly anisotropic diffusion tensor; in the isotropic rank-2 model its spectral density has Lorentzian frequency dependence proportional to `tau_c/[1 + (omega*tau_c)^2]`.

After creating and basing the system, the function evaluates `relaxation(spin_system)` and prints its full matrix. It defines no initial or detection operator and produces no observable or spectrum, so this example demonstrates the relaxation model and superoperator, not a simulated experimental signal.
