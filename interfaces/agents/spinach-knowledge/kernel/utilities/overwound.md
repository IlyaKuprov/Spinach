# kernel/utilities/overwound.m

- Signature: `overwound(rho,spc_dim,spn_dim)`

## Purpose

Checks if a Fokker-Planck state vector has any spatial frequencies that its spatial grid is dangerously close to misrepresenting due to insufficient point count. Syntax: overwound(rho,spc_dim,spn_dim)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Parameters / inputs

- rho -Fokker-Planck state vector or
- a bookshelf stack thereof
- spc_dim -spatial dimensions of the Fokker-
- Planck problem, [X Y Z]
- spn_dim -spin dimension of the Fokker-
- Planck problem
- Output:
- figures and diagnostic messages to the console
- Note: if you are running spatial dynamics, such as diffusion and
- and flow, with finite difference derivative operators, the
- spatial grid point count shuld be set to several times the
- minimum Nyquist value.

## Implementation structure

- Checks if a Fokker-Planck state vector has any spatial frequencies
- that its spatial grid is dangerously close to misrepresenting due
- to insufficient point count. Syntax:
- overwound(rho,spc_dim,spn_dim)
- rho -Fokker-Planck state vector or
- a bookshelf stack thereof
- spc_dim -spatial dimensions of the Fokker-
- Planck problem, [X Y Z]
- spn_dim -spin dimension of the Fokker-
- Planck problem
- Output:
- figures and diagnostic messages to the console
