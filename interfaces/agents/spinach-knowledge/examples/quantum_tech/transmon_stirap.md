# examples/quantum_tech/transmon_stirap.m

- Signature: `transmon_stirap()`

## Purpose

Ensemble-GRAPE optimization of a single-transmon STIRAP transfer using Gaussian pulses and a penalty on excess power. Calculation time: minutes. The source identifies the model with https://doi.org/10.1038/ncomms10628.

## Model and parameters

The rotating-frame T3 transmon has a 10 MHz ladder detuning and -20 MHz anharmonicity. Gaussian pulses drive the 0–1 and 1–2 transitions over 300 points spanning -150 to 150 ns, with 45 ns width and a -90 ns pulse-pair delay. Their nominal normalized amplitudes are set by 43.4 and 38.2 MHz, and the pulse-power ensemble is 30, 35, 40, 45, and 50 MHz. The transfer is from BL1 to BL3.

## Optimization

The two transition controls are optimized by GRAPE with the SNSA penalty (weight 1.0), Goodwin method, and a 50-iteration limit. The Gaussian pulse pair is the initial guess; the source cites https://doi.org/10.1038/ncomms10628.
