# etc/textbook/rlx_sbm.m

## Signature

`[r1,r2]=rlx_sbm(B0,nucleus,dist,a_iso,e_spin,g_eff,t1e,t2e,tau_r)`

## Purpose

Calculates nuclear longitudinal and transverse relaxation rates caused by a paramagnetic centre, using the Solomon–Bloembergen–Morgan (SBM) model. Dipolar and isotropic-hyperfine contact mechanisms are reported separately.

## Model and calculation

The dipolar terms depend on the electron–nucleus distance, the effective electron spin and g-factor, and correlation times combining rotational motion with the electron's longitudinal or transverse relaxation. The contact terms depend on the isotropic hyperfine coupling and electron relaxation times. The spectral-density convention is `J(omega,tau)=tau/(1+omega^2*tau^2)`.

## Inputs

- `B0`: positive magnetic field in tesla.
- `nucleus`: nuclear isotope label, for example `'1H'` or `'13C'`.
- `dist`: electron–nucleus distance in angstroms.
- `a_iso`: isotropic hyperfine coupling in radians per second.
- `e_spin`: effective electron spin quantum number.
- `g_eff`: effective electron g-factor.
- `t1e`, `t2e`: electron longitudinal and transverse relaxation times, respectively, in seconds.
- `tau_r`: rotational correlation time in seconds.

## Outputs

Each output is a two-element vector `[dipolar, contact]`:

- `r1`: longitudinal relaxation rates in Hz.
- `r2`: transverse relaxation rates in Hz.
