# examples/optimal_control/bloch_siegert/yusuke_optimal_vs_cw_demo.m

- Signature: `yusuke_optimal_vs_cw_demo()`
- Status: historical documentation. The corresponding `.m` file is absent from the current checkout; this description is based on `b4f03f29^:examples/optimal_control/bloch_siegert/yusuke_optimal_vs_cw_demo.m` and does not establish current behavior.

## Purpose

The historical example optimises a Bloch–Siegert-aware, phase-modulated low-power waveform against a constant-phase X-pulse baseline. It is a reduced single-spin identity-cycle surrogate intended to preserve magnetisation over offset and RF-amplitude (B1) variation. The header describes it as a control-side companion to `yusuke_14n_broadening_demo.m`; it is explicitly not a full QJF/MAS quadrupolar simulation.

## Historical setup and numerics

- Random seed: `rng(1)`.
- Field: `sys.magnet=18.8` T, corresponding to 800 MHz for (^1)H.
- Model: one (^ {13})C spin, zero isotropic Zeeman coupling, and the full spherical-tensor Liouville basis (`formalism='sphten-liouv'`, `approximation='none'`). The identity-cycle states are normalised (S_x), (S_y), and (S_z).
- Header context: 70 kHz MAS, (t_psim10) μs, and (
u_{14N}) in the 15–23 kHz range. These are motivating 14N-decoupling conditions, not a simulated quadrupolar 14N system.
- RF and pulse: nominal RF 20 kHz; 10 equal elements of 10 μs each (100 μs total). The constant-phase reference is a 4π X pulse. The phase-only optimiser uses L-BFGS, up to 40 iterations, with Bloch–Siegert corrections enabled (`control.bsiegert=true()`); amplitudes are held fixed.
- Optimisation training grid: seven offsets from −12 to +12 kHz and B1 scales [0.95, 1.00, 1.05]. Evaluation grid: 61 offsets from −20 to +20 kHz and nine B1 scales from 0.90 to 1.10.

The objective evaluates preservation of the three Cartesian basis states; the code also evaluates offset and B1 profiles for the optimised and constant-phase waveforms and plots phase/profile comparisons. The source contains no prose conclusion or fixed numerical result; this page therefore makes no claim that either waveform wins over the evaluation grid.

## Source and citation status

The historical header says the numerical regime is inspired by “the 14N decoupling papers” but gives no bibliographic citation. No paper citation or DOI is supplied by the source. The file's author line is `aditya.dev@weizmann.ac.il`.
