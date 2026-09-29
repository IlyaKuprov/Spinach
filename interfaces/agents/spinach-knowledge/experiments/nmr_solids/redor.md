# experiments/nmr_solids/redor.m

MATLAB source: [experiments/nmr_solids/redor.m](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/nmr_solids/redor.m)

REDOR (rotational-echo double resonance) returns the full echo, dephased echo, and their difference as a function of rotor-cycle count. The source describes ideal hard `pi` pulses on the dephasing channel and is intended for a `singlerot` context. Its examples include `13C{15N}`; the cited references are [10.1016/0022-2364(89)90280-1](https://doi.org/10.1016/0022-2364(89)90280-1) and [10.1006/jmre.2000.2128](https://doi.org/10.1006/jmre.2000.2128).

## Inputs and source-defined sequence

`redor_curve=redor(spin_system,parameters,H,R,K)` forms `L=H+1i*R+1i*K` and analytically decouples any nuclei listed in `parameters.decouple` (default `{}`). `parameters.spins` orders the observed and dephasing spins; the pulse operators are constructed on the second spin. The transverse operators are extended by `speye(parameters.spc_dim)`, where `spc_dim` comes from the context. This function receives, rather than builds, the Hamiltonian and spatial/rotor context.

Required fields are `spc_dim`, `spins`, `ncycles`, `rate`, `rho0`, and `coil`. `rate` is the MAS rate in Hz, so the code uses `1/abs(rate)` seconds per rotor period and half that duration for each half-period. `ncycles` selects the requested rotor-cycle samples. `decouple` is an optional cell array of isotope strings and must not include either REDOR working spin. `pulse_phase` is an optional nonempty real row vector; it defaults to `[0 pi/2]` and is used as phase angles in `cos(phase)*Ix+sin(phase)*Iy`.

The reference echo propagates under `L` for a full rotor period per cycle. The dephased echo propagates for two half-periods per cycle, with a hard `pi` pulse after each half-period; successive pulses use `pulse_phase` cyclically. Both echoes are detected as `coil'*rho`, including at zero cycles. The code retains samples from zero through `max(ncycles)` and then selects the requested counts.

## Output and limits

The result has three rows and one column per requested cycle count: row 1 is `S0`, row 2 is `S`, and row 3 is `S0-S`. The source requires `sphten-liouv` formalism and compatible matrix inputs for `H`, `R`, and `K`. The function specifies rotor-synchronised ideal pulses, not finite RF pulse shapes; orientation and the underlying MAS Hamiltonian are supplied by the calling context rather than swept here.

https://spindynamics.org/wiki/index.php?title=redor.m
