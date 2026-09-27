# examples/nmr_stochastic/snmr_gb1.m

- Signature: `snmr_gb1()`

## Purpose

Runs a Primas-style stochastic NMR trajectory for GB1 protein. The source warns that the calculation requires a terabyte of RAM and an NVIDIA A100 GPU and estimates hours. Author: Ilya Kuprov (ilya.kuprov@weizmann.ac.il).

## Model and calculation

The protein is imported from 2N9K.pdb and 2N9K.bmrb with all selections retained, then simulated at 18.79 T. The basis uses IK-1, scalar-coupling connectivity, interaction level 4, and proximity level 3; interaction and proximity cutoffs are 2.0 and 4.0. Redfield relaxation uses kite retention, IME equilibrium, correlation time 5e-9 s, and temperature 298 K.

The script builds the Hamiltonian and relaxation superoperator, initializes isotropic thermal equilibrium, and defines x/y control and observable operators for 1H, 13C, and 15N. It generates six independent Gaussian control-noise tracks with sigma=100 Hz, timestep 1e-5 s, and 1e6 steps. GPU arrays are used for the operators and state. At each step it records six expectation values and advances the state with the time-dependent generator; checkpoints are saved every 1000 steps to gb1_workspace.mat.

## Output

The example reports the trajectory rate and plots each control track alongside its corresponding observable trajectory.
