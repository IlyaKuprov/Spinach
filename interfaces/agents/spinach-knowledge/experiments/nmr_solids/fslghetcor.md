# experiments/nmr_solids/fslghetcor.m

- Signature: `fid=fslghetcor(spin_system,parameters,H,R,K)`

## Purpose

Heteronuclear correlation MAS NMR experiment with frequency-switched Lee-Goldburg (FSLG) homonuclear decoupling. Further details: https://doi.org/10.1103/PhysRev.140.A1261, https://doi.org/10.1016/0009-2614(89)87166-0, and https://doi.org/10.1006/jmre.1996.1089 (Figure 1).

## Parameters / inputs

- `spin_system` — spin system; the implementation requires `sphten-liouv` formalism.
- `parameters.spins` — two working spin isotopes, e.g. `{'1H','13C'}`; the first is the high-gamma channel.
- `parameters.hi_pwr` — amplitude of high-power pulses on the high-gamma channel, Hz.
- `parameters.cp_pwr` — pulse amplitudes on the two channels during the cross-polarisation (CP) contact time, Hz.
- `parameters.cp_dur` — CP contact duration, s.
- `parameters.offset` — transmitter offsets on the two channels, Hz.
- `parameters.nblocks` — number of FSLG blocks per indirect-dimension point.
- `parameters.spc_dim` — Fokker-Planck spatial dimension.
- `parameters.rho0` — initial state.
- `parameters.coil` — detection state.
- `parameters.sweep` — `[F1 F2]` sweep widths, Hz. The F1 element is unused because the F1 dwell time is set by the FSLG block duration; it may be `NaN`.
- `parameters.npoints` — numbers of points in F1 and F2.
- `H` — Hamiltonian superoperator, received from the context function.
- `R` — relaxation superoperator, received from the context function.
- `K` — kinetics superoperator, received from the context function.

## Outputs

- `fid.sin`, `fid.cos` — sine and cosine components of the States quadrature.

## Implementation

The sequence generates separate cosine and sine F1 trajectories using alternating FSLG evolution blocks, returns the first channel from the magic angle, applies the CP contact, decouples that channel during F2 acquisition, and separates the two quadrature components. Control operators are extended across `parameters.spc_dim`; GPU execution is used for FSLG generators when enabled.

Source: <https://spindynamics.org/wiki/index.php?title=fslghetcor.m>