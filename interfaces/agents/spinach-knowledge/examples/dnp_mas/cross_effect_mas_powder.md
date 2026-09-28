# examples/dnp_mas/cross_effect_mas_powder.m

- Signature: `cross_effect_mas_powder()`

## Purpose

Calculates the steady-state DNP enhancement for a powder under MAS using Spinach's `masdnp` routine, following [Mentink-Vigier et al.](http://dx.doi.org/10.1016/j.jmr.2015.07.001). The source notes that Spinach uses different rotation conventions from the paper and gives a runtime of minutes.

## Physical and numerical setup

The model comprises two electron spins and one proton at 9.394 T, with anisotropic electron g tensors, electron-electron and electron-proton couplings, and Nottingham relaxation. It specifies electron and nuclear longitudinal/transverse rates, 100 K temperature, Di Bari equilibrium, and secular relaxation retention; the basis is the complete spherical-tensor Liouville basis.

The `masdnp` calculation uses the electron spin, a 12.5 kHz rotor rate, the specified MAS axis, rank 800, microwave power parameter `2*pi*0.85e6`, microwave frequency `-263.366e9`, and a 1.0 s microwave time. Powder orientations use `rep_2ang_100pts_sph`. The proton coil state is set to its `Lz` state. The routine's returned enhancement is displayed as the steady-state DNP enhancement.
