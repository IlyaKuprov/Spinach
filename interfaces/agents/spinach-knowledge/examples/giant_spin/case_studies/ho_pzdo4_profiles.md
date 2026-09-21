# examples/giant_spin/case_studies/ho_pzdo4_profiles.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/giant_spin/case_studies/ho_pzdo4_profiles.m`
- Signature: `ho_pzdo4_profiles()`
- Total lines: 117

## Purpose

Pulsed-field magnetisation of the Ho(pzdo)4 metal-organic framework, a J=8 giant spin with a crystal field to twelfth spherical rank, under four magnetic field profiles: linear sweep, piecewise linear sweep, monotone cubic spline through a measured 65 T short pulse, and a sinusoidal field at the clock transition frequency. Spin-phonon relaxation is the generalised Lindblad dissipator of Saito and Miyashita with a super-Ohmic phonon bath. Reproduces Figure 2 of https://arxiv.org/abs/2609.16352 with the crystal field parameters, g-factor, temperatures, spectral density, sweep profiles, and stair widths of that paper. As in the paper, the first three profiles are propagated for 1 ms (the first millisecond of the measured 10 ms pulse), the sinusoid for 140 ps. Calculation time: minutes.

## Physical / mathematical content

- Giant-spin examples. The effective model treats lanthanides or high-spin centres using crystal-field / Stevens-operator Hamiltonians, Zeeman splitting, and magnetisation dynamics.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- Pulsed-field magnetisation of the Ho(pzdo)4 metal-organic framework,
- a J=8 giant spin with a crystal field to twelfth spherical rank, under
- four magnetic field profiles: linear sweep, piecewise linear sweep,
- monotone cubic spline through a measured 65 T short pulse, and a sinu-
- soidal field at the clock transition frequency. Spin-phonon relaxation
- is the generalised Lindblad dissipator of Saito and Miyashita with a
- super-Ohmic phonon bath. Reproduces Figure 2 of
- with the crystal field parameters, g-factor, temperatures, spectral
- density, sweep profiles, and stair widths of that paper. As in the
- paper, the first three profiles are propagated for 1 ms (the first
- millisecond of the measured 10 ms pulse), the sinusoid for 140 ps.
- Calculation time: minutes

## Internal Spinach / MATLAB structure cues

- Called routines in the main body: `ho_pzdo4_params()`, `stev2sph()`, `icm2hz()`, `stevens()`, `double()`, `pchip()`, `kfigure()`, `scale_figure()`, `create()`, `basis()`, `crystal()`, `subplot()`, `plot()`, `kylabel()`.
