# examples/nmr_stochastic/snmr_strychnine.m

- Signature: `snmr_strychnine()`

## Purpose

A Primas-style stochastic NMR experiment on strychnine. The calculation requires an NVidia Titan V GPU at a minimum. Calculation time: minutes

## Physical / mathematical content

- Stochastic NMR examples. These scripts model random processes, trajectories, or stochastic Liouville dynamics and connect fluctuating Hamiltonians or transport processes to ensemble-averaged observables.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Implementation structure

- A Primas-style stochastic NMR experiment on strychnine. The calculation
- requires an NVidia Titan V GPU at a minimum.
- Calculation time: minutes
- Read the spin system properties
- Magnet field
- Algorithmic options
- Basis set
- Relaxation theory parameters
- Use GPU arithmetic
- sys.enable={'gpu'};
- Spinach housekeeping
- Get the Hamiltonian
