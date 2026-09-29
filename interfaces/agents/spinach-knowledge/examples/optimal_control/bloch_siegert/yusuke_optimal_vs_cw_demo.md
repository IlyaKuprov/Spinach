# examples/optimal_control/bloch_siegert/yusuke_optimal_vs_cw_demo.m

- Signature: `yusuke_optimal_vs_cw_demo()`

## Aim and model boundary

The script optimises a phase-modulated, low-power identity cycle and compares it with a constant-phase X-pulse having the same RF amplitude and total duration. The intended reduced-model question is whether the cycle preserves three Cartesian magnetisation states across offset and B1 variation. The source calls this a control-side companion to `yusuke_14n_broadening_demo.m`, but explicitly frames it as a surrogate—not a full quadrupolar-(^{14}mathrm N), QJF, or MAS simulation. Its header gives 70 kHz MAS, approximately 10 μs pulse elements, and a 15–23 kHz (^{14}mathrm N) frequency range as motivation only.

## Historical setup

The script fixes `rng(1)`, sets `sys.magnet=18.8` T (800 MHz for (^{1}mathrm H)), and creates a one-spin (^{13}mathrm C) model with zero isotropic Zeeman coupling. It uses `formalism='sphten-liouv'` and `approximation='none'`. Normalised (S_x,S_y,S_z) states serve as the three initial and target operators; (L_x,L_y) are the two RF controls and (L_z) the offset operator.

The pulse has 10 equal 10 μs elements (100 μs total) at nominal 20 kHz RF. Amplitudes remain fixed while element phases are varied. The constant-phase baseline uses zero phase for every element; the header characterises it as a 4π X pulse. The optimiser is configured for L-BFGS with at most 40 iterations, uses Bloch–Siegert corrections, and starts from a random phase vector. The training set comprises seven offsets from −12 to +12 kHz and B1 scales ([0.95,1.00,1.05]).

## Evaluation and interpretation

Both waveforms are evaluated over 61 offsets from −20 to +20 kHz and nine B1 scaling factors from 0.90 to 1.10. The script produces fidelity maps and offset/B1 profiles, including mean profiles, for visual comparison. These grids describe the script's evaluation procedure; they are not reported experimental results. The source contains no prose conclusion or fixed numerical result, so it does not establish that either waveform outperforms the other across the evaluation range.

The historical header mentions “the 14N decoupling papers” but names no paper or DOI. No bibliographic citation should be inferred from the filename or motivation.
