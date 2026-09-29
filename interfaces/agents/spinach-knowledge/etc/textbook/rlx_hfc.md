# etc/textbook/rlx_hfc.m

## Use

[r1,r2,rx]=rlx_hfc(B0,HFC,spins,tau_c) computes Redfield relaxation rates for two spins coupled by a hyperfine tensor during isotropic tumbling in a liquid. Use the documented pair as one electron and one nuclear spin; the implementation expects labels understood by Spinach's spin helper.

## Inputs

- B0: real scalar magnetic field in tesla.
- HFC: real 3-by-3 hyperfine coupling tensor in rad/s. It need not be symmetric.
- spins: two character isotope/spin labels in a cell array, in the order used for outputs; e.g. {'E','15N'}.
- tau_c: positive real scalar rotational correlation time in seconds.

## Calculation and outputs

The code obtains rank-1 and rank-2 Blicharski invariants of HFC, spin-square factors and the two Zeeman frequencies, then combines the invariants with isotropic-tumbling spectral densities at zero, single-spin, sum, and difference frequencies. The rotational diffusion parameter passed to the spectral-density calls is 1/(6*tau_c). The rank contributions enter longitudinal and transverse rates separately; rx combines the difference-frequency rank-1 term with rank-2 sum- and difference-frequency terms.

- r1: two longitudinal rates in input-spin order, Hz.
- r2: two transverse rates in input-spin order, Hz.
- rx: longitudinal cross-relaxation rate, Hz.

## Scope and source

The routine is for the stated isotropic-tumbling Redfield model, not a general motion model. Input checks enforce real scalar B0, a real 3-by-3 tensor, two character labels, and positive real scalar tau_c; they do not explicitly check that one label is an electron, so supply the intended electron–nucleus pair. Source: [implementation](https://github.com/IlyaKuprov/Spinach/blob/main/etc/textbook/rlx_hfc.m) · [Spin Dynamics Wiki](https://spindynamics.org/wiki/index.php?title=rlx_hfc.m).
