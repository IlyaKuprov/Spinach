# etc/textbook/rlx_sbm.m

## Use

[r1,r2]=rlx_sbm(B0,nucleus,dist,a_iso,e_spin,g_eff,t1e,t2e,tau_r) estimates nuclear relaxation from a paramagnetic centre with the Solomon–Bloembergen–Morgan model. Each returned row vector has separate dipolar and isotropic-contact contributions.

## Inputs

- B0: positive finite magnetic field in tesla.
- nucleus: nuclear isotope label understood by Spinach's spin helper, e.g. '1H' or '13C'.
- dist: positive finite electron–nucleus distance in Å.
- a_iso: finite real isotropic hyperfine coupling in rad/s.
- e_spin: positive integer or half-integer effective electron spin quantum number.
- g_eff: positive effective electron g-factor.
- t1e, t2e: positive finite electron longitudinal and transverse relaxation times, respectively, in seconds.
- tau_r: positive finite rotational correlation time in seconds.

These are scalar inputs; the numerical parameters are checked as real and finite where applicable.

## Calculation and outputs

The dipolar calculation forms electron and nuclear Zeeman frequencies, combines tau_r with each electron relaxation time to obtain separate effective correlation times, and uses the dipolar prefactor proportional to dist^-6. The spectral-density convention documented in the source is J(omega,tau) = tau/(1 + omega^2*tau^2). The contact terms depend on a_iso, e_spin, t1e, t2e, and the electron–nuclear frequency difference.

- r1 = [r1_dd, r1_sc]: longitudinal dipolar and contact rates, Hz.
- r2 = [r2_dd, r2_sc]: transverse dipolar and contact rates, Hz.

## Scope and source

The model uses an isotropic hyperfine contact term and the specified distance-based dipolar term; it does not take an anisotropic hyperfine tensor. The isotope string must also resolve through spin. Source: [implementation](https://github.com/IlyaKuprov/Spinach/blob/main/etc/textbook/rlx_sbm.m) · [Spin Dynamics Wiki](https://spindynamics.org/wiki/index.php?title=rlx_sbm.m).
