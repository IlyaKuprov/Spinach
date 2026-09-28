# kernel/utilities/offsetof.m

- Signature: `offs=offsetof(spin_system,idx)`

## Purpose

Returns the isotropic Zeeman offset of the specified spin from the pure magnetogyric-ratio frequency in the current magnet.

## Parameters

- `spin_system` — spin system containing the Zeeman tensors and base frequencies.
- `idx` — index of the spin in the isotope array. Use `idxof()` to find the index from a text label. Must be a positive integer no greater than the number of particles in `spin_system.comp.isotopes`.

## Output

- `offs` — offset from the pure magnetogyric-ratio frequency at the current field, in Hz.

## Calculation

After checking `idx`, the function takes the spin’s Zeeman tensor from `spin_system.inter.zeeman.matrix{idx}`, subtracts `eye(3)*spin_system.inter.basefrqs(idx)`, and converts its isotropic part to Hz: `offs=-trace(offs)/(3*2*pi)`.

Source: <https://spindynamics.org/wiki/index.php?title=offsetof.m>