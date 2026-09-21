# examples/giant_spin/case_studies/ho_pzdo4_profiles.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/giant_spin/case_studies/ho_pzdo4_profiles.m`
- Signature: `ho_pzdo4_profiles()`
- Total lines: 147

## Purpose

Pulsed-field magnetisation of the Ho(pzdo)4 metal-organic framework, a J=8 giant spin with a crystal field to twelfth spherical rank, under four magnetic field profiles: linear sweep, piecewise linear sweep, monotone cubic spline through a measured 65 T short pulse, and a sinusoidal field at the clock transition frequency. Spin-phonon relaxation is the generalised Lindblad dissipator of Saito and Miyashita with a super-Ohmic phonon bath. Reproduces Fig 2 of https://arxiv.org/abs/2609.16352 with the crystal field parameters, g-factor, temperatures, spectral density, sweep profiles, stair widths, panel layout, and axis limits of that paper. As in the paper, the first three profiles are propagated for 1 ms (the first millisecond of the measured 10 ms pulse), the sinusoid for ten periods of 46.8 ps. Calculation time: minutes.

## Physical / mathematical content

- One `E17` spin (J=8) with g=1.24 and the crystal field of `ho_pzdo4_params.m`, single crystal with the crystal field frame aligned with the laboratory frame.
- The four `parameters.field_prof` handles are 10 T/ms linear, `interp1` piecewise linear through the paper's turning points (1 us, 0.1 T), (10 us, 1 T), (100 us, 5 T), and (1 ms, 10 T), `pchip` through the 76 measured points of the paper's 65 T pulse, and a 0.1 T sinusoid at the angular frequency of the 0.71 cm^-1 clock gap (period 46.8 ps); the first three run at 2 K, the sinusoid at 1 mK where only the clock doublet is populated and the magnetisation oscillates non-periodically between about -4.7 and +5.4 mu_B.
- The spin-phonon coupling operator has unit elements between adjacent m_J states; super-Ohmic bath (`phonon_alpha=2`) with lambda=10 cm^-1 and I0=1e-14 ps/rad; the observable is the moment -1.24*J_z in Bohr magnetons.

## Numerical / algorithmic content

- `zeeman-hilb`, `crystal` context with the `labframe` assumption and `parameters.needs={'zeeman_op'}`, `sys.magnet=1`. Stair widths are 10 ns for the millisecond profiles (10^5 stairs, one record every 100) and 0.1 ps for the sinusoid (4680 stairs, one record every 4), the paper's values.
- The spin system is re-created for every panel because the bath temperature changes, and J_z (`operator(spin_system,'Lz','E17')`), the coupling mask, and the coil are built after `basis` inside the loop.
- Panels (a) to (d) as in the paper: magnetisation in red on the left axis (0 to 6 mu_B, -6 to 6 for the sinusoid), field in black on the right axis (0 to 10 T, 0 to 13 T for the spline, -0.25 to 0.25 T for the sinusoid), time in picoseconds with the paper's ticks, and the temperature and time step written in the lower right corner.
