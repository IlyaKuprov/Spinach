# experiments/nmr_solids/fslghetcor.m

Source: [canonical MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/nmr_solids/fslghetcor.m) · [Spin Dynamics Wiki](https://spindynamics.org/wiki/index.php?title=fslghetcor.m)

Heteronuclear-correlation MAS NMR with frequency-switched Lee–Goldburg (FSLG) homonuclear decoupling and a cross-polarisation (CP) contact. The source cites [10.1103/PhysRev.140.A1261](https://doi.org/10.1103/PhysRev.140.A1261), [10.1016/0009-2614(89)87166-0](https://doi.org/10.1016/0009-2614(89)87166-0), and [10.1006/jmre.1996.1089](https://doi.org/10.1006/jmre.1996.1089) (Figure 1).

## Inputs

Signature: `fid=fslghetcor(spin_system,parameters,H,R,K)`. It requires the `sphten-liouv` formalism; `H`, `R`, and `K` must be same-sized numeric matrices.

- `parameters.spins` has exactly two isotope strings present in the spin system. The source example is the common `1H`/`13C` pair; the first entry supplies the FSLG controls and the second the second-channel control operator.
- `hi_pwr` is a positive scalar high-power-channel amplitude in Hz; `cp_pwr` is a two-element vector of positive CP-channel amplitudes in Hz; `cp_dur` is a positive contact duration in seconds.
- `offset` is a two-element real vector in Hz. `nblocks` and `spc_dim` are positive integers. `rho0` and `coil` are required initial and detection states.
- `sweep` is a two-element real vector in Hz; its F1 element may be NaN because F1 timing comes from the FSLG blocks. `npoints` is a two-element vector of positive integers for F1 and F2.

## Sequence outline

The code extends the `Lx`, `Ly`, and `Lz` controls over the Fokker–Planck dimension, and forms `L=H+1i*R+1i*K`. It starts cosine and sine trajectories with the high-power pulse for `1/(4*hi_pwr)` plus the magic-angle duration `acos(1/sqrt(3))/(2*pi*hi_pwr)`. Alternating FSLG generators use the first offset and `hi_pwr/sqrt(2)`; each block applies paired evolutions of `sqrt(2/3)/hi_pwr`. The implementation optionally moves those generators to a GPU when GPU support is enabled.

It then applies the paired CP generator `L-2*pi*cp_pwr(1)*Hy+2*pi*cp_pwr(2)*Cx` for `cp_dur`, decouples the first listed spin, and acquires F2 with `coil` at dwell `1/sweep(2)`. Output fields `fid.cos` and `fid.sin` are the two States-quadrature components; each axis uses its corresponding `npoints` count. The source does not define additional protein-specific transfer or REDOR steps.