# examples/dnp_mas/solid_effect_mas_enlev.m

- MATLAB implementation: [examples/dnp_mas/solid_effect_mas_enlev.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/dnp_mas/solid_effect_mas_enlev.m)

## Purpose

This no-argument example, run as `solid_effect_mas_enlev()`, builds an ESR rotor stack for an electron–`^1H` pair and plots the energy levels as the rotor phase changes. It follows Fred Mentink-Vigier's MAS DNP treatment, with the source noting that Spinach's rotation conventions differ ([paper](https://doi.org/10.1016/j.jmr.2015.07.001)). The source header estimates milliseconds.

## Model and stack

The field is assigned `9.403`. The electron g-tensor eigenvalues are `[2.00614 2.00194 2.00988]` with Euler angles `pi*[253.6 105.1 123.8]/180`; the proton Zeeman eigenvalues and Euler angles are both zero. Coordinates are `[0 0 0]` and `[0 0 3.00]`. The basis is full Zeeman-Hilbert (`formalism='zeeman-hilb'`, `approximation='none'`).

The rotor-stack axis is `[sqrt(2/3) 0 sqrt(1/3)]`, with empty rotor frames, orientation `[0 0 0]`, spins `{'E','1H'}`, magnet MAS frame, zero offsets, and `max_rank=200`. The code calls `rotor_stack(spin_system,parameters,'esr')`.

## Output and scope

For each stack Hamiltonian, the example sorts the real eigenvalues of `eig(H{n})`, then plots them against rotor phase over `[0,2*pi)`. The axes are rotor phase in radians and level energy in rad/s. This is a single-orientation energy-level plot; the code does not perform powder averaging or propagate a microwave-driven, relaxing state.
