# examples/giant_spin/case_studies/ho_pzdo4_profiles.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/giant_spin/case_studies/ho_pzdo4_profiles.m`
- Signature: `ho_pzdo4_profiles()`
- Total lines: 114

## Purpose

Pulsed-field magnetisation of the Ho(pzdo)4 metal-organic framework, a J=8 giant spin with a crystal field to twelfth spherical rank, under four magnetic field profiles: linear sweep, piecewise linear sweep, monotone cubic spline through a measured 65 T short pulse, and a sinusoidal field at the clock transition frequency. Spin-phonon relaxation is the generalised Lindblad dissipator of Saito and Miyashita with a super-Ohmic phonon bath. Reproduces Fig 2 of https://arxiv.org/abs/2609.16352 with the crystal field parameters, g-factor, temperatures, spectral density, sweep profiles, and stair widths of that paper. As in the paper, the first three profiles are propagated for 1 ms (the first millisecond of the measured 10 ms pulse), the sinusoid for 140 ps. Calculation time: minutes.

## Physical / mathematical content

- One `E17` spin (J=8) with g=1.24 and the crystal field of `ho_pzdo4_params.m`, single crystal with the crystal field frame aligned with the laboratory frame.
- The four `parameters.field_prof` handles are 10 T/ms linear, `interp1` piecewise linear through five points to 10 T at 1 ms, `pchip` through the 76 measured points of the paper's 65 T pulse, and a 0.1 T sinusoid at the angular frequency of the clock gap; the first three run at 2 K, the sinusoid at 1 mK where only the clock doublet is populated.
- The spin-phonon coupling operator has unit elements between adjacent m_J states; super-Ohmic bath (`phonon_alpha=2`) with lambda=10 cm^-1 and I0=1e-14 ps/rad; the observable is the moment -1.24*J_z in Bohr magnetons.

## Numerical / algorithmic content

- `zeeman-hilb`, `crystal` context with the `labframe` assumption and `parameters.needs={'zeeman_op'}`, `sys.magnet=1`. Stair widths are 10 ns for the millisecond profiles (10^5 stairs, one record every 100) and 1 fs for the sinusoid (140538 stairs, one record every 59); the plots use microseconds and picoseconds accordingly.
- The spin system is re-created for every panel because the bath temperature changes, and J_z (`operator(spin_system,'Lz','E17')`), the coupling mask, and the coil are built after `basis` inside the loop.
