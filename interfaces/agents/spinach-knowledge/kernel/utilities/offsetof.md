# kernel/utilities/offsetof.m

## Purpose

Returns the isotropic Zeeman offset of a specified spin from the pure magnetogyric ratio frequency in the current magnet.

## Behaviour

- Syntax: `offs=offsetof(spin_system,idx)`.
- Validates the spin index via an internal consistency check (`grumble`): the index must be a positive integer scalar (numeric, real, integral, at least 1), otherwise an error `'idx must be a positive integer.'` is raised.
- Errors with `'idx exceeds the number of particles in the system.'` if the index exceeds `numel(spin_system.comp.isotopes)`.
- Retrieves the Zeeman tensor from `spin_system.inter.zeeman.matrix{idx}`.
- Subtracts the magnet frequency term: `offs - eye(3)*spin_system.inter.basefrqs(idx)`.
- Extracts the isotropic part as `trace(offs)/3` and converts to Hz via `offs = -offs/(2*pi)`.

## Inputs and outputs

**Inputs**

- `spin_system` — spin system object.
- `idx` — index of the spin in the `sys.isotopes` array; use `idxof()` to find the index by the text label.

**Outputs**

- `offs` — offset from the pure magnetogyric ratio frequency at the current field, in Hz.

## References

- Source: [kernel/utilities/offsetof.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/offsetof.m)
- Spinach Wiki: [offsetof.m](https://spindynamics.org/wiki/index.php?title=offsetof.m)
