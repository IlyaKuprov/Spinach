# examples/giant_spin/case_studies/ho_twelfth_rank/dimer_exchange_types.m

- MATLAB implementation: [examples/giant_spin/case_studies/ho_twelfth_rank/dimer_exchange_types.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/giant_spin/case_studies/ho_twelfth_rank/dimer_exchange_types.m)

- Signature: `dimer_exchange_types()`
- Source: `examples/giant_spin/case_studies/ho_twelfth_rank/dimer_exchange_types.m`

## Purpose

Compares the pulsed-field magnetisation of an electron-spin dimer under four exchange tensors with the corresponding thermal-equilibrium magnetisation. The example reproduces Fig. 4 of [arXiv:2609.16352](https://arxiv.org/abs/2609.16352), including its colours, line styles, and plot limits. The setup is a single crystal at 0.2 K swept to 1 T at 10 T/ms; the reported calculation time is minutes.

## Physical interpretation

The two spins are S=1/2 with g=2. Each run uses one tensor: isotropic 0.2 cm^-1, Jxx-only, Jzz-only, or the antisymmetric tensor shown in the source. The paper writes the interaction as H=-2 S1·J·S2; Spinach uses S1·A·S2 with A in hertz, so the conversion is A=-2 icm2hz(J).

The isotropic and Jzz-only interactions commute with the Zeeman term and do not coherently mix Zeeman levels. In those cases population transfer is due to the spin-phonon dissipator, which couples product states whose total S_z differs by one. The Jxx-only and antisymmetric tensors mix levels, allowing the swept magnetisation to approach equilibrium with a lag. The existing entry reports the 1 T sweep/equilibrium values as approximately 0.006/2.0 μ_B (isotropic), 0.003/2.0 μ_B (Jzz), 0.61/1.98 μ_B (Jxx), and 1.61/1.90 μ_B (antisymmetric); these are retained from that entry, not newly measured here.

The bath is super-Ohmic (phonon_alpha=2), using the paper's lambda^2 I0 with lambda=10 cm^-1 and I0=1e-10 ps/rad, converted to rad/s units. The spin-phonon coupling mask selects transitions between adjacent total-S_z states, and the plotted observable is the total moment -2 S_z in μ_B.

## Calculation and output

The simulation uses the zeeman-hilb formalism and a crystal context with the labframe assumption. It computes a thermal-equilibrium curve at each field for comparison with the pulsed-field trajectory; the equilibrium values are retained in each result as `answers{n}.obs_eq`. The sweep uses 10 ns steps for 10^4 steps and records every tenth step. The plot overlays equilibrium (solid) and sweep (dashed) curves in red, blue, green, and orange for tensors 1–4, with field limits 0–1 T and magnetisation limits 0–2 μ_B. The legend distinguishes “J_n^dimer Equilibrium” and “J_n^dimer QME”.
