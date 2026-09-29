# experiments/bruker/hmqcetgp.m

- Source: [experiments/bruker/hmqcetgp.m](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/bruker/hmqcetgp.m)
- Signature: `fid=hmqcetgp(spin_system,parameters,H,R,K)`

## Purpose and sequence model

Simulate echo/antiecho gradient-selected HMQC, following the Bruker `hmqcetgp` pulse program and a standard HMQC sequence. The source comment cites the standard-sequence paper at [DOI 10.1016/0022-2364(83)90241-X](https://doi.org/10.1016/0022-2364(83)90241-X); this is a sequence citation, not evidence that a particular simulated spectrum or experimental result was validated.

The routine begins with longitudinal magnetisation on the F2 spin, applies a `pi/2` pulse on F2, evolves for `abs(1/(2*J))` seconds (with `J` supplied in Hz), and applies pulses on F1. During the indirect F1 evolution it builds a trajectory, applies configured midpoint refocusing pulses, and separates the `+1` and `-1` F1 coherence-order pathways. Subsequent pulse/coherence-selection operations form the echo and antiecho branches; each is detected on the F2 `L+` state during direct-dimension acquisition. Gradient selection is represented by coherence-order selection, not explicit gradient pulses.

The sequence uses the caller-supplied `L = H + 1i*R + 1i*K`; the function does not build the chemical spin system or assign its Zeeman, scalar-coupling, relaxation, or kinetic parameters. Any such model is encoded in the input matrices. The implementation requires the `sphten-liouv` formalism.

## Parameters / inputs

- `parameters.sweep`: two sweep widths `[F1 F2]`, in Hz.
- `parameters.npoints`: point counts `[F1 F2]`.
- `parameters.spins`: F1 and F2 isotope strings, for example `{'13C','1H'}`.
- `parameters.J`: working scalar coupling in Hz, used to set the transfer delay `abs(1/(2*J))`.
- `parameters.decouple_f1`: nuclei receiving midpoint 180-degree refocusing pulses in F1, e.g. `{'1H'}`.
- `parameters.decouple_f2`: nuclei to decouple in F2, e.g. `{'15N','13C'}`.
- `H`, `R`, and `K`: Hamiltonian, relaxation, and kinetics matrices from the context function.

The routine uses dwell times `1/sweep(1)` and `1/sweep(2)` for F1 and F2, respectively.

## Outputs

- `fid.pos` and `fid.neg`: echo and antiecho signal components.

The source recommends isotope dilution for natural-abundance simulations; see [dilute.m](https://spindynamics.org/wiki/index.php?title=dilute.m). This routine's observable is the F2 detected HMQC signal, not a singlet-yield observable.

## References and links

- [Standard HMQC sequence cited in the source](https://doi.org/10.1016/0022-2364(83)90241-X)
- [Spinach documentation for hmqcetgp.m](https://spindynamics.org/wiki/index.php?title=hmqcetgp.m)
