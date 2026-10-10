# examples/relaxation_theory/aniso_diff_test_1.m

[Source](https://github.com/IlyaKuprov/Spinach/blob/main/examples/relaxation_theory/aniso_diff_test_1.m)

This example constructs a Redfield relaxation superoperator for a heteronuclear 1H–13C pair in anisotropic rotational diffusion. It does not simulate an acquisition or compare a predicted signal with experimental data: its result is the full matrix R, printed to the console.

The field is assigned as 2*pi*950.33e6/spin('1H'). Both chemical-shielding tensors are set to zero, so the fluctuating interaction specified by the geometry is the dipolar coupling: the spins are separated by 1.13 Å, with the vector rotated by euler2dcm(1,2,3). The static scalar coupling is 145.0 Hz. The rotational-diffusion tensor eigenvalues are [2.16e8 2.35e8 7.45e8]; the corresponding correlation-time parameters supplied to Spinach are 1./(6*D). The file does not tabulate a spectral-density curve.

Relaxation is requested with the Redfield model and zero equilibrium. It retains the lab-frame relaxation matrix and uses the untruncated sphten-liouv basis. No particular cross-correlation pair is explicitly selected or separately reported. The calculation calls relaxation(spin_system) and displays full(R); it reports neither R1/R2 projections nor an experimental observable.
