# kernel/summaries/summary_pbc.m

- Signature: `summary_pbc(spin_system,header)`

## Behaviour

Prints each periodic-boundary vector from `spin_system.inter.pbc` as X, Y and Z columns. Each component is signed and formatted to three decimal places; the source assigns no coordinate units. The supplied `header` precedes the table.

The function returns nothing; every line is routed through `report(spin_system,...)`. The local guard requires a structure `spin_system` and character-array `header`; it does not validate individual vector lengths.

## References

- [MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/summaries/summary_pbc.m)
- [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=summary_pbc.m)
