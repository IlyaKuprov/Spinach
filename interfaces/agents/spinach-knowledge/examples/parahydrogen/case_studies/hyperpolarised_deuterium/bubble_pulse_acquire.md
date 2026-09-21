# examples/parahydrogen/case_studies/hyperpolarised_deuterium/bubble_pulse_acquire.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/parahydrogen/case_studies/hyperpolarised_deuterium/bubble_pulse_acquire.m`
- Signature: `bubble_pulse_acquire()`
- Total lines: 120

## Purpose

Simulated PNL (partially negative line) spectrum of ortho-deuterium in the presence of a parahydrogena- tion catalyst. Bubbling is followed by a 45-degree pulse. Paper link to follow in due course.

## Physical / mathematical content

- Parahydrogen examples. The physical motif is highly non-Boltzmann singlet order imported from para-H2 and converted into observable nuclear magnetisation through hydrogenation, exchange, or catalytic transfer processes.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- Simulated PNL (partially negative line) spectrum of
- ortho-deuterium in the presence of a parahydrogena-
- tion catalyst. Bubbling is followed by a 45-degree
- pulse. Paper link to follow in due course.
- Spin system
- Experimental chemical shifts
- Experimental J-couplings
- NQI tensors from a DFT calculation
- Cartesian coordinates, DFT calculation
- Kinetics
- Magnet field
- Simulation formalsim

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `relaxan()`, `deut_pair()`, `unit_state()`, `assume()`, `hamiltonian()`, `relaxation()`, `kinetics()`, `pumped_state()`, `magpump()`, `evolution()`, `state()`, `operator()`, `liquid()`, `apodisation()`.
