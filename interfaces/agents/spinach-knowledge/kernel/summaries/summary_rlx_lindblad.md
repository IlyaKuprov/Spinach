# kernel/summaries/summary_rlx_lindblad.m

- Signature: `summary_rlx_lindblad(spin_system,header)`

## Behaviour

Prints one Lindblad relaxation-rate row per spin, in `spin_system.comp.nspins` order. Each row contains its index, isotope, `spin_system.rlx.lind_r1_rates(n)`, `spin_system.rlx.lind_r2_rates(n)`, and spin label; scalar rates use signed scientific notation with five decimal places. No rate-unit conversion occurs in this routine.

The supplied header precedes the table and all output is sent through `report`; there is no return value. The local guard requires a structure `spin_system` and character-array `header`.

## References

- [MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/summaries/summary_rlx_lindblad.m)
- [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=summary_rlx_lindblad.m)
