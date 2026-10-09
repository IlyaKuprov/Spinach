# examples/dnp_sol/crosspol_powder_static_1.m

- MATLAB implementation: [examples/dnp_sol/crosspol_powder_static_1.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/dnp_sol/crosspol_powder_static_1.m)

## Purpose

This no-argument example models a static-powder, doubly rotating-frame electron–`15N` cross-polarisation contact experiment. It computes and plots the nitrogen `S_x` expectation value over the contact pulse. The source estimates the calculation time as seconds.

## Model and sequence

The system contains `15N` and an electron (`sys.isotopes={'15N','E'}`), with `sys.magnet=9.394`, Zeeman scalars 0 and 2.0023193043622, coordinates `[0,0,0]` and `[0,0,10.05]`, and `temperature=298`. No coordinate units are given in the script. It uses the full `sphten-liouv` basis (`approximation='none'`) and does not assign a relaxation model in this example.

The sequence has 100 intervals, each `1e-5` seconds, so the plotted time axis runs from zero to 1 ms with 101 samples. Both rows of `irr_powers` contain 100 values of `5e4`. The supplied irradiation operators are electron `Ly` and nitrogen `Lx`; the excitation operators are electron `Lx` and nitrogen `Ly`. The detection state is nitrogen `Lx`, `spins={'15N'}`, and the powder grid is `rep_2ang_6400pts_sph`. The script requests `needs={'iso_eq'}` and comments that this is “Good enough here”; that qualification belongs to this example's chosen equilibrium treatment.

## Calculation and output

Spinach's `powder` driver calls `@cp_contact_hard` in NMR mode: `powder(spin_system,@cp_contact_hard,parameters,'nmr')`. The returned local variable `fid` is plotted as `real(fid)` against cumulative contact time in seconds; the vertical axis is labelled as the nitrogen `S_x` expectation value. The function itself declares no output argument, so its result is the generated figure rather than a returned FID. Running the example requires Spinach and the `cp_contact_hard` callback.
