# examples/kinetics/relayed_hyperpol.m

- Signature: `relayed_hyperpol()`

## Purpose

Relayed NOE from hyperpolarized water to an ALA–GLY dipeptide, generating Figure S7 from Christopher Pötzl: https://doi.org/10.1016/j.jmr.2024.107727.

## Physical / mathematical content

The 30-proton model contains ten molecular protons and 20 water protons. The water spins have no coordinates, preventing direct cross-relaxation in this model; intermolecular exchange connects them to the first four labile protons at 20 s⁻¹. Redfield relaxation with empirical water R1 and R2 rates of 0.1 Hz, an IK-1 basis retaining up to three-spin molecular orders, and the specified thermal-equilibrium state are used. The initial state is changed to 100% water polarisation.

## Numerical / algorithmic content

The combined Hamiltonian, relaxation, and exchange-kinetics Liouvillian is propagated for 128 steps of 0.125 s. Detection tracks the Z magnetisation of the three aliphatic protons and the α proton, and plots their time courses.

## Implementation structure

- Sets B₀ = 16.4 T and 30 protons; water shifts are 4.5 ppm.
- Uses exchange flux rate 20 between protons 1–4 and water spins 11–20.
- Builds `L=H+1i*R+1i*K`, prepares full water polarisation, and calls `evolution` in multichannel mode.
- Plots the CH₃ and Hα magnetisation in arbitrary units over 16 s.
