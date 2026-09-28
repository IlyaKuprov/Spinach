# kernel/spin.m

- Signature: `[gamma,multiplicity]=spin(name)`

## Purpose

Returns the magnetogyric ratio and multiplicity for a named isotope or supported particle, including spin-zero species.

## Physical / mathematical content

- `gamma` is the magnetogyric ratio in rad/(s*Tesla).
- `multiplicity` is the number of energy levels or population levels.

## Numerical / algorithmic content

The function looks up isotope data by name. It also handles the ghost-spin case `G` (`gamma=0`, multiplicity 1), neutron `N`, and muon `M`. The parameterized forms `E#`, `C#`, `V#`, and `T#` represent, respectively, a high-spin electron, electromagnetic cavity mode, phonon mode, and transmon; `#` specifies the multiplicity or number of levels. `E#` requires at least two levels; the cavity, phonon, and transmon forms require at least three. Their magnetogyric ratio is zero for cavity, phonon, and transmon modes.

## Parameters / inputs

- `name` — isotope or supported-particle name, such as `'15N'` or `'195Pt'`; special forms are listed above.

## Outputs

- `gamma` — magnetogyric ratio in rad/(s*Tesla); zero for cavities, phonons, and transmons.
- `multiplicity` — number of energy or population levels.

## Notes

The source warns that entries without a stated source were sourced from Google and should be double-checked before production calculations. Some known isotopes have no data in the current NMR literature; other unrecognized names produce an unknown-isotope error.
