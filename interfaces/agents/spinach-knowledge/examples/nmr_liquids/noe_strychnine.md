# examples/nmr_liquids/noe_strychnine.m

- Signature: `noe_strychnine()`

## Purpose

An inversion-recovery NOE-effect spectrum calculation for the strychnine proton system. Spin 9 is inverted, allowed to relax for 500 ms, and the resulting difference from equilibrium is measured with pulse-acquire. The source estimates a calculation time of minutes.

## Spin system and relaxation model

The example obtains the system from `strychnine({'1H'})`, sets `sys.magnet=14.1`, and disables Krylov propagation. It uses the spherical-tensor Liouville formalism, the IK-2 approximation, scalar-coupling connectivity, proximity level 3, and a proximity cut-off of 4.0. Relaxation is Redfield; the equilibrium convention is `dibari`, retained terms are `kite`, temperature is set to 298, and `tau_c={200e-12}`.

## Preparation, acquisition, and processing

The code constructs the Redfield relaxation superoperator and thermal equilibrium state. It forms the inverted state as `rho_eq-2*Lz9*(Lz9'*rho_eq)/norm(Lz9)^2`, propagates it under `1i*R` for 0.5, then subtracts equilibrium to leave the perturbation.

Pulse-acquire observes `1H`: the initial state is the difference operator, the coil is `L+`, and a `Ly` pulse of angle `pi/2` is applied. Decoupling is empty. The source sets offset 2800, sweep 6500, 8192 points, zero-filling to 65536, ppm axis units, and an inverted axis. It runs `liquid(...,@hp_acquire,...,'nmr')`, applies exponential apodisation with parameter 6, Fourier-transforms and shifts the result, then plots the real spectrum.

## Source

[examples/nmr_liquids/noe_strychnine.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_liquids/noe_strychnine.m)
