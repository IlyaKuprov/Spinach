# examples/relaxation_theory/dd_csa_xcorr_2.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/relaxation_theory/dd_csa_xcorr_2.m`
- Signature: `dd_csa_xcorr_2()`
- Total lines: 103

## Purpose

DD-CSA cross-correlation -a reproduction of Fig 5a from the paper by Grace and Kumar (http://dx.doi.org/10.1006/jmra.1995.1151). Calculation time: seconds

## Physical / mathematical content

- Relaxation-theory examples. The mathematical backbone is Bloch-Redfield-Wangsness or stochastic Liouville theory, spectral densities, cross-correlation terms, motional models, and extraction of longitudinal/transverse decay behaviour from superoperators.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.
- Chemical-shift anisotropy is present: shielding is treated as a second-rank tensor whose orientation relative to the field or rotor axis modulates line shapes and transfer dynamics.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- DD-CSA cross-correlation -a reproduction of Fig 5a from the paper
- by Grace and Kumar (http://dx.doi.org/10.1006/jmra.1995.1151).
- Calculation time: seconds
- Read the spin system parameters (vacuum DFT calculation)
- Set up the calculation
- Proximity cut-off
- Run Spinach housekeeping
- Set simulation parameters
- Set the assumptions to high-field NMR
- Get the Hamiltonian superoperator
- Add Redfield superoperator,
- Apply the offset

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `g2spinach()`, `gparse()`, `create()`, `basis()`, `assume()`, `hamiltonian()`, `relaxation()`, `frqoffset()`, `operator()`, `state()`, `sweep2ticks()`, `equilibrium()`, `kfigure()`, `step()`, `evolution()`, `fftshift()`.
