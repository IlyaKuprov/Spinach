# examples/quantum_tech/transmon_transfer.m

Source: [examples/quantum_tech/transmon_transfer.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/quantum_tech/transmon_transfer.m)

## What it models

GRAPE optimisation of coherence transfer from transmon 1 to transmon 2 in a coupled pair of Duffing modes. This is a simulated control problem with power and frequency-offset ensembles, not an experimental transfer result; the source estimates calculation time in minutes.

## Drift model and states

The source sets the field to zero and declares T3 and T5 modes with rotating-frame frequencies 100e6 and -200e6 Hz, anharmonicities -10e6 and -20e6 Hz, and an inter-mode exchange parameter of 50e6 Hz. The Zeeman-Hilbert basis is unapproximated. The cavity drift includes the two Duffing ladders and their flip-flop exchange coupling.

The initial state is the sum of C and A coherence operators on transmon 1 with transmon 2 in BL1; the target places the same C/A coherence on transmon 2 while transmon 1 is in BL1. Both are Frobenius-normalised, symmetrised, and the target is rescaled by the Sorensen bound, as specified by the source comment. C and A denote the raising and lowering operators used to form the coherence components.

## Control ensemble and optimisation

Two quadrature controls address the respective transmons. Each has five offset samples spanning -10e6 to 10e6 Hz, and the source sets the power ensemble to 2*pi*[40,45,50,55,60]*1e6*5 rad/s. The 200 pulse slices are 0.25 ns each. The random two-channel initial pulse is optimised with GRAPE using the SNSA penalty (weight 1.0), rbfgs method, and 200-iteration limit; control, spectrogram, and robustness plots are enabled. The source contains no relaxation model, rotational diffusion, correlation spectrum, cross-correlations, secular approximation, or reported final fidelity. It therefore supports no relaxation-rate benchmark or experimental-agreement claim.
