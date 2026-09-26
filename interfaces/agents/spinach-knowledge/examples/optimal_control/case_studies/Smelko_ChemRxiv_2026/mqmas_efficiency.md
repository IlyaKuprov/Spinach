# examples/optimal_control/case_studies/Smelko_ChemRxiv_2026/mqmas_efficiency.m

- Signature: `mqmas_efficiency()`

## Purpose

Efficiency of the z-filtered 27Al MQMAS pulse sequence with hard pulses and with the optimal control pulses produced by the other examples in this folder. Reproduces, using Spinach, the sequence efficiency calculation from https://doi.org/10.26434/chemrxiv.15008427 The nucleus has the quadrupolar coupling of aluminium acetylacetonate (CQ=3.2 MHz, eta=0.16) and is spun at 12.5 kHz in a 400 MHz magnet; the quadrupolar interaction is taken to second order in the rotating frame. The sequence is: excitation pulse, +MQ/-MQ coherence filter, conversion pulse, population filter, central transition selective pulse, and detection of the central transition single-quantum coherence. The efficiency is the modulus of the detected element of the density matrix, normalised to the initial state as in the paper. The powder average runs over 400 crystallite orientations at 32 initial rotor phases each. Hard pulses have the durations optimised in the paper. Optimal control waveforms are read from the files written by mq_excitation.m, mq_conversion.m, and ct_selective.m examples for both coherence orders; the files supplied in this folder were produced by those examples on 128 cores and only need rerunning if the optimisations are changed. With these waveforms, the 3QMAS efficiency is 0.071 with hard pulses and 0.419 with optimal control pulses, and the 5QMAS efficiency is 0.0109 and 0.258, respectively: signal enhancement factors of 5.9 and 23.8, against 5.7 and 25 simulated in the paper.

## Physical / mathematical content

- Optimal-control examples. These scripts formulate pulse design as a nonlinear optimisation problem over waveform samples or basis coefficients. The core mathematical objects are fidelities, gradients, Hessians or Hessian approximations, ensemble robustness objectives, and constrained search over RF amplitude/phase trajectories.
- Quadrupolar physics is relevant: nuclei with spin > 1/2 interact with the electric field gradient tensor, introducing second-rank anisotropy, asymmetry, and overtone or MQ phenomena.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Implementation structure

- Efficiency of the z-filtered 27Al MQMAS pulse sequence with hard
- pulses and with the optimal control pulses produced by the other
- examples in this folder. Reproduces, using Spinach, the sequence
- efficiency calculation from
- The nucleus has the quadrupolar coupling of aluminium acetylace-
- tonate (CQ=3.2 MHz, eta=0.16) and is spun at 12.5 kHz in a 400
- MHz magnet; the quadrupolar interaction is taken to second order
- in the rotating frame. The sequence is: excitation pulse, +MQ/-MQ
- coherence filter, conversion pulse, population filter, central
- transition selective pulse, and detection of the central transi-
- tion single-quantum coherence. The efficiency is the modulus of
- the detected element of the density matrix, normalised to the
