# experiments/bruker/hsqcetgpsi.m

- Signature: `fid=hsqcetgpsi(spin_system,parameters,H,R,K)`

## Purpose

Sensitivity-improved, echo/antiecho gradient-selected HSQC based on the Bruker `hsqcetgpsi` pulse program and standard HSQC. Gradient selection is represented by coherence-order selection statements. The routine returns the two acquisition pathways as `fid.pos` and `fid.neg`.

## Implementation

The routine builds `L=H+1i*R+1i*K`, sets the evolution delay from the working scalar coupling `parameters.J`, and prepares the proton-detected HSQC pathways using the specified F1 and F2 spins. It selects the pathway coherences analytically and acquires the echo and antiecho signals. For natural-abundance simulations, the source recommends isotope-dilution functionality (see `dilute.m`).

## Parameters / inputs

- `parameters.sweep`: [F1 F2] sweep widths, Hz.
- `parameters.npoints`: [F1 F2] numbers of points.
- `parameters.spins`: {F1 F2} nuclei, e.g. {'13C','1H'}.
- `parameters.decouple_f2`: nuclei to decouple in F2, e.g. {'15N','13C'}.
- `parameters.decouple_f1`: nuclei receiving midpoint 180-degree refocusing pulses in F1, e.g. {'1H','15N'}; must not include the F1 active isotope.
- `parameters.J`: working scalar coupling, Hz.
- `parameters.trim_angle`: proton trim-pulse angle, rad.
- `parameters.si_time`: sensitivity-improvement delay, s.
- `H`: Hamiltonian matrix supplied by the context function.
- `R`: relaxation superoperator supplied by the context function.
- `K`: kinetics superoperator supplied by the context function.

## Output

- `fid.pos`, `fid.neg`: detected echo and antiecho signal components.

## References

The source cites these standard-HSQC references:

- [10.1016/0009-2614(80)80041-8](https://doi.org/10.1016/0009-2614(80)80041-8)
- [10.1002/cmr.a.10095](https://doi.org/10.1002/cmr.a.10095)
- [10.1016/0022-2364(91)90036-S](https://doi.org/10.1016/0022-2364(91)90036-S)
- [10.1021/ja00052a088](https://doi.org/10.1021/ja00052a088)
- [10.1007/BF00175254](https://doi.org/10.1007/BF00175254)

[Source page](https://spindynamics.org/wiki/index.php?title=hsqcetgpsi.m)
