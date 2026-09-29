# examples/dnp_mas/cross_effect_mas_enlev.m

- MATLAB implementation: [examples/dnp_mas/cross_effect_mas_enlev.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/dnp_mas/cross_effect_mas_enlev.m)

- Function: `cross_effect_mas_enlev()`

## Purpose

Builds the three-spin cross-effect MAS model and plots its sorted Hamiltonian energy levels over one rotor period. The example follows the MAS DNP model associated with [Mentink-Vigier et al.](https://doi.org/10.1016/j.jmr.2015.07.001); its source notes that Spinach rotation conventions differ from the paper.

## Model and rotor stack

The spins are `{'E','E','1H'}`, with `sys.magnet=9.394`. The electron Zeeman principal values are `[2.0094 2.0060 2.0017]` for both electrons; the first tensor uses zero Euler angles and the second uses `pi*[107 108 124]/180`. The proton Zeeman values are zero. Electron–electron coupling values are `[23.0e6 -11.5e6 -11.5e6]` with Euler angles `pi*[0 135 0]/180`; the first electron–proton coupling is `[1.5e6 -0.75e6 -0.75e6]` with zero Euler angles.

The complete Zeeman-Hilbert basis (`formalism='zeeman-hilb'`, `approximation='none'`) is used. `rotor_stack(spin_system,parameters,'labframe')` uses axis `[sqrt(2/3) 0 sqrt(1/3)]`, rate `12.5e3`, orientation `pi*[320 141 80]/180`, magnetic MAS frame, selected spins `{'E','1H'}`, zero offsets, empty rotor frames, and `max_rank=200`.

## Result and scope

For every stack Hamiltonian, the code sorts the real eigenvalues and divides by `2*pi*1e9`; the plot labels energy in GHz and time in μs. It shows levels 1–2, 3–6, and 7–8 in three panels, with the eigensolves parallelised using `parfor`. The source estimates milliseconds of calculation time. This is an energy-level plot, not a density-operator trajectory or a calculation of DNP enhancement; it includes no relaxation or microwave-drive term.
