# examples/nmr_stochastic/snmr_gb1.m

- Signature: `snmr_gb1()`

## Purpose

A Primas-style stochastic NMR experiment on GB1 protein. The calculation requires a terabyte of RAM and NVidia A100 GPU. Calculation time: hours

## Physical / mathematical content

- Stochastic NMR examples. These scripts model random processes, trajectories, or stochastic Liouville dynamics and connect fluctuating Hamiltonians or transport processes to ensemble-averaged observables.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Implementation structure

- A Primas-style stochastic NMR experiment on GB1 protein. The
- calculation requires a terabyte of RAM and NVidia A100 GPU.
- Calculation time: hours
- Protein data import
- Magnet field
- Tolerances
- Basis set
- Relaxation theory
- Use GPU arithmetic
- sys.enable={'gpu'};
- Spinach housekeeping
- Get the Hamiltonian
