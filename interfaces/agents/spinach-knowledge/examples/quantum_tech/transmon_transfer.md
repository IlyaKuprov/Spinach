# examples/quantum_tech/transmon_transfer.m

- Signature: `transmon_transfer()`

## Purpose

Coherence transfer from transmon 1 to transmon 2 in a coupled two-transmon Duffing model. GRAPE optimization accounts for distributions of control powers and transmon offsets and penalizes excess power. Calculation time: minutes.

## Model and parameters

The T3 and T5 modes have rotating-frame frequencies 100 MHz and -200 MHz, anharmonicities -10 MHz and -20 MHz, and a 50 MHz flip-flop exchange coupling. The initial coherence is on transmon 1 and the target coherence is on transmon 2.

## Optimization

The two control channels use five offset samples each, spanning -10 to 10 MHz per transmon. The pulse-power levels are `2*pi*[40,45,50,55,60] MHz*5`; the 200 slices are 0.25 ns each. GRAPE uses the SNSA penalty (weight 1.0), rbfgs method, and a 200-iteration limit, starting from a random two-channel pulse. The target is normalized using the Sørensen bound.
