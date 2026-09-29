# kernel/utilities/chemshifts.m

## Purpose

Returns the chemical shifts of every spin in the system relative to the carrier frequency in the current magnet. Syntax:

```
[cs_ppm,cs_hz]=chemshifts(spin_system)
```

## Behaviour

- The function first calls an internal consistency check (`grumble`) that errors with `'the spin system object does not contain the required information.'` if the `spin_system` object lacks either the `comp` or the `inter` field.
- Outputs `cs_ppm` and `cs_hz` are preallocated as `spin_system.comp.nspins`-by-1 zero vectors.
- For each spin `n` from 1 to `spin_system.comp.nspins`:
  - The isotropic Zeeman frequency is obtained as `iso=trace(spin_system.inter.zeeman.matrix{n})/3`.
  - The carrier frequency is subtracted: `iso=iso-spin_system.inter.basefrqs(n)`.
  - The chemical shift in ppm is computed as `cs_ppm(n)=1e6*iso/spin_system.inter.basefrqs(n)`.
  - The chemical shift in Hz is computed as `cs_hz(n)=-iso/(2*pi)`.

## Inputs and outputs

**Inputs**

- `spin_system` — spin system descriptor object; must contain the `comp` and `inter` fields.

**Outputs**

- `cs_ppm` — chemical shifts in ppm, one per spin.
- `cs_hz` — chemical shifts in Hz, one per spin.

## References

- Source: [chemshifts.m on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/chemshifts.m)
- Spin Dynamics Wiki: <https://spindynamics.org/wiki/index.php?title=chemshifts.m>
