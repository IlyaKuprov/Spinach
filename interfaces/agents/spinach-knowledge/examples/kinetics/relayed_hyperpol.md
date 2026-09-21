# examples/kinetics/relayed_hyperpol.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/kinetics/relayed_hyperpol.m`
- Signature: `relayed_hyperpol()`
- Total lines: 101

## Purpose

Relayed NOE from hyperpolarized water to ALA-GLY dipeptide, generating Figure S7 from Christopher Pötzl

## Physical / mathematical content

- Chemical-kinetics examples. The files couple spin dynamics to exchange, pumping, or nonlinear reaction networks represented by kinetic generators in Liouville space.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

## Implementation structure

- Relayed NOE from hyperpolarized water to ALA-GLY dipeptide,
- generating Figure S7 from
- Christopher Pötzl
- Simulation timing parameters
- Magnet field
- 30 protons in the system
- Cartesian coordinates of pertinent protons
- 20 water protons exist, but have no coordinates to
- prevent direct cross-relaxation from happening
- Chemical shifts, all water at 4.5 ppm
- Relaxation theories
- Empirical relaxation at 0.1 Hz for water

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `repelem()`, `num2cell()`, `create()`, `basis()`, `assume()`, `hamiltonian()`, `relaxation()`, `kinetics()`, `equilibrium()`, `state()`, `evolution()`, `figure()`, `scale_figure()`, `kxlabel()`, `kylabel()`, `klegend()`.
