# examples/nmr_stochastic/snmr_gb1.m

- Signature: `snmr_gb1()`
- Source: [`examples/nmr_stochastic/snmr_gb1.m`](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_stochastic/snmr_gb1.m)

## Experiment and model

This is a Primas-style stochastic NMR simulation of GB1, not a measured trajectory. The source credits Ilya Kuprov. It imports molecular input from `2N9K.pdb` and `2N9K.bmrb`, selects molecule 1 and all atoms, and deletes shifts lacking data. The system is set to 18.79 T. The sphten-liouv basis uses IK-1 with scalar-coupling connectivity, interaction level 4 and proximity level 3; the source also sets interaction and proximity cutoffs to 2.0 and 4.0. Relaxation is Redfield, with the kite retention option, IME equilibrium, correlation time 5 ns and temperature 298 K.

## Stochastic drive and signal

The script constructs the NMR Hamiltonian and relaxation superoperator, then evolves the thermal-equilibrium state under the base generator H0 + iR and six independently generated Gaussian control tracks: x/y controls for 1H, 13C and 15N. Each track has sigma = 2π × 100 rad/s (100 Hz) and is sampled at dt = 10 μs for 1,000,000 steps, a nominal 10 s trajectory. At every step it records the six corresponding Lx/Ly expectation values and advances the state with the time-dependent generator. Checkpoints are saved every 1,000 steps to `gb1_workspace.mat`; the script also saves the Hamiltonian, relaxation operator and equilibrium workspace.

The output is a 6-by-2 plot pairing each control in Hz with its calculated expectation-value trajectory. There is no pulse, gradient, chirp or field-drop schedule in this example. GPU arrays are used in the trajectory section, while the explicit `sys.enable={'gpu'}` line is commented out. The source warns of terabyte-scale memory, an NVIDIA A100 and hours of calculation; those are source comments, not independently verified resource measurements.
