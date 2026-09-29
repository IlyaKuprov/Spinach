# kernel/summaries/summary_symmetry.m

- Signature: `summary_symmetry(spin_system,header)`

## Behaviour

Reports a two-column table headed Group and Spins. Row `n` prints `spin_system.comp.sym_group{n}` and its associated spin-index list `spin_system.comp.sym_spins{n}`; it does not calculate the symmetry groups. The supplied `header` precedes the table.

The function returns nothing and sends its lines through `report`. The local guard requires a structure `spin_system` and character-array `header`.

## References

- [MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/summaries/summary_symmetry.m)
- [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=summary_symmetry.m)
